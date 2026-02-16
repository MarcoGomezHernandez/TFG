#include "genetic.hpp"

/*
 * Genetic algorithm configuration constants (private)
 */
namespace GeneticPrivateConfig
{
    // Population and generations
    static constexpr size_t POPULATION_SIZE = 50;
    static constexpr size_t NUM_GENERATIONS = 200;
    static constexpr size_t NUM_ELITES = 2;

    // Observation time is use_time / this divisor
    static constexpr double OBSERVATION_TIME_DIVISOR = 5.0;

    // Crossover and mutation probabilities
    static constexpr double CROSSOVER_PROBABILITY = 0.8;
    static constexpr double MUTATION_PROBABILITY = 0.1;

    // Random initialization ranges [min, max] for each mutable parameter
    static constexpr double ESYN_MIN = -2.0;
    static constexpr double ESYN_MAX = 2.0;
    static constexpr double SFAST_MIN = 0.1;
    static constexpr double SFAST_MAX = 5.0;
    static constexpr double VFAST_MIN = -2.0;
    static constexpr double VFAST_MAX = 2.0;
    static constexpr double VSLOW_MIN = -2.0;
    static constexpr double VSLOW_MAX = 2.0;
    static constexpr double K1_MIN = 0.001;
    static constexpr double K1_MAX = 1.0;
    static constexpr double K2_MIN = 0.001;
    static constexpr double K2_MAX = 1.0;
    static constexpr double SSLOW_MIN = 0.1;
    static constexpr double SSLOW_MAX = 5.0;

    // Mutation perturbation factors (multiplied by N(0,1)) per parameter
    static constexpr double ESYN_MUT_FACTOR = 0.2;
    static constexpr double SFAST_MUT_FACTOR = 0.3;
    static constexpr double VFAST_MUT_FACTOR = 0.2;
    static constexpr double VSLOW_MUT_FACTOR = 0.2;
    static constexpr double K1_MUT_FACTOR = 0.05;
    static constexpr double K2_MUT_FACTOR = 0.05;
    static constexpr double SSLOW_MUT_FACTOR = 0.3;
}

/*
 * Comparator for descending fitness order (higher fitness first)
 */
static bool fitness_descending(const Individual &a, const Individual &b)
{
    return a.fitness > b.fitness;
}

/*
 * Generate initial population with random individuals
 */
template <size_t POP_SIZE>
inline std::array<Individual, POP_SIZE> initialize_population(std::mt19937 &rng)
{
    std::uniform_real_distribution<double> dist_esyn(GeneticPrivateConfig::ESYN_MIN, GeneticPrivateConfig::ESYN_MAX);
    std::uniform_real_distribution<double> dist_sfast(GeneticPrivateConfig::SFAST_MIN, GeneticPrivateConfig::SFAST_MAX);
    std::uniform_real_distribution<double> dist_vfast(GeneticPrivateConfig::VFAST_MIN, GeneticPrivateConfig::VFAST_MAX);
    std::uniform_real_distribution<double> dist_vslow(GeneticPrivateConfig::VSLOW_MIN, GeneticPrivateConfig::VSLOW_MAX);
    std::uniform_real_distribution<double> dist_k1(GeneticPrivateConfig::K1_MIN, GeneticPrivateConfig::K1_MAX);
    std::uniform_real_distribution<double> dist_k2(GeneticPrivateConfig::K2_MIN, GeneticPrivateConfig::K2_MAX);
    std::uniform_real_distribution<double> dist_sslow(GeneticPrivateConfig::SSLOW_MIN, GeneticPrivateConfig::SSLOW_MAX);

    std::array<Individual, POP_SIZE> population;

    for (Individual &ind : population)
    {
        ChemicalSynapsisParams &p = ind.params;
        p.Esyn = dist_esyn(rng);
        p.sfast = dist_sfast(rng);
        p.Vfast = dist_vfast(rng);
        p.Vslow = dist_vslow(rng);
        p.k1 = dist_k1(rng);
        p.k2 = dist_k2(rng);
        p.sslow = dist_sslow(rng);
    }

    return population;
}

/*
 * Arithmetic crossover: operate directly on ChemicalSynapsisParams (child = mean of two parents)
 */
inline void crossover(const ChemicalSynapsisParams &a, const ChemicalSynapsisParams &b, ChemicalSynapsisParams &result)
{
    result.Esyn = (a.Esyn + b.Esyn) / 2.0;
    result.sfast = (a.sfast + b.sfast) / 2.0;
    result.Vfast = (a.Vfast + b.Vfast) / 2.0;
    result.Vslow = (a.Vslow + b.Vslow) / 2.0;
    result.k1 = (a.k1 + b.k1) / 2.0;
    result.k2 = (a.k2 + b.k2) / 2.0;
    result.sslow = (a.sslow + b.sslow) / 2.0;
}

/*
 * Mutate: recibe directamente los parámetros de la sinapsis (params) y aplica
 * perturbaciones normales por cada campo con su probabilidad independiente.
 */
inline void mutate(ChemicalSynapsisParams &p, std::mt19937 &rng, std::normal_distribution<double> &ndist, std::uniform_real_distribution<double> &prob_dist)
{
    constexpr double MUTATION_PROBABILITY = GeneticPrivateConfig::MUTATION_PROBABILITY;
    if (prob_dist(rng) < MUTATION_PROBABILITY)
        p.Esyn += ndist(rng) * GeneticPrivateConfig::ESYN_MUT_FACTOR;
    if (prob_dist(rng) < MUTATION_PROBABILITY)
        p.sfast += ndist(rng) * GeneticPrivateConfig::SFAST_MUT_FACTOR;
    if (prob_dist(rng) < MUTATION_PROBABILITY)
        p.Vfast += ndist(rng) * GeneticPrivateConfig::VFAST_MUT_FACTOR;
    if (prob_dist(rng) < MUTATION_PROBABILITY)
        p.Vslow += ndist(rng) * GeneticPrivateConfig::VSLOW_MUT_FACTOR;
    if (prob_dist(rng) < MUTATION_PROBABILITY)
        p.k1 += ndist(rng) * GeneticPrivateConfig::K1_MUT_FACTOR;
    if (prob_dist(rng) < MUTATION_PROBABILITY)
        p.k2 += ndist(rng) * GeneticPrivateConfig::K2_MUT_FACTOR;
    if (prob_dist(rng) < MUTATION_PROBABILITY)
        p.sslow += ndist(rng) * GeneticPrivateConfig::SSLOW_MUT_FACTOR;
}

/*
 * Roulette wheel selection: returns reference to selected individual
 * roulette_dist is precomputed with range [0, fitness_sum]
 * ignore_ind: pointer to individual to skip during selection (nullptr to ignore none)
 */
template <size_t N>
inline const Individual &roulette_select_one(const std::array<Individual, N> &population,
                                             std::uniform_real_distribution<double> &roulette_dist,
                                             std::mt19937 &rng)
{
    double pick = roulette_dist(rng);
    double cumulative = 0.0;
    for (const Individual &ind : population)
    {
        cumulative += ind.fitness;
        if (cumulative >= pick)
            return ind;
    }
    // Fallback sencillo: devolver el último individuo
    return population.back();
}

/*
 * Roulette wheel selection excluding a specific individual
 * Uses the same distribution but skips the ignored individual during accumulation
 */
template <size_t N>
inline const Individual &roulette_select_second_no_rep(const std::array<Individual, N> &population,
                                                       const Individual &ignore_ind,
                                                       std::uniform_real_distribution<double> &roulette_dist,
                                                       std::mt19937 &rng)
{
    double pick = roulette_dist(rng);
    double cumulative = 0.0;
    for (const Individual &ind : population)
    {
        if (&ind == &ignore_ind)
            continue;
        cumulative += ind.fitness;
        if (cumulative >= pick)
            return ind;
    }

    // Fallback: return first individual that is not the ignored one
    for (const Individual &ind : population)
    {
        if (&ind != &ignore_ind)
            return ind;
    }
    return population.back();
}

/*
 * Main genetic algorithm function template.
 * NeuronType: the neuron wrapper type
 * Integrator: the numeric integrator type
 * CreateFuncType: callable returning NeuronType
 * ResetStateFuncType: callable taking NeuronType& and resetting state
 * GetVFuncType: callable taking const NeuronType& and returning double (voltage)
 */
template <typename Integrator, typename NeuronType,
          CreateFunc<NeuronType> CreateFuncType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
Individual genetic(const std::string &csv_path,
                   size_t column_index,
                   double csv_step,
                   double start_time,
                   double use_time,
                   NumericIntegrator integrator,
                   NeuronModel model,
                   bool search_phase,
                   bool check_drift,
                   CreateFuncType create_neuron,
                   ResetStateFuncType reset_state_neur,
                   GetVFuncType get_v_neur)
{
    // Calculate observation time
    const double observation_time = use_time / GeneticPrivateConfig::OBSERVATION_TIME_DIVISOR;

    // --- Step 1: Scale the signal ---
    ScaledSignalResult scaled_result = scale_signal(
        csv_path, column_index, csv_step, start_time, use_time,
        observation_time, integrator, model, check_drift);

    if (!scaled_result.success)
    {
        throw std::runtime_error("Signal scaling failed.");
    }

    // --- Step 2: Create neuron and synapsis instances ---
    NeuronType model_neur = create_neuron();
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;
    ChemicalSynapsisType synapsis;
    // Set fixed synapsis parameters once
    synapsis.set(ChemicalSynapsisType::gfast, GeneticPublicConfig::GFAST_FIXED);
    synapsis.set(ChemicalSynapsisType::gslow, GeneticPublicConfig::GSLOW_FIXED);

    // --- Step 3: Allocate buffers ---
    const size_t signal_size = scaled_result.signal.size();
    std::vector<double> model_signal_buffer;
    model_signal_buffer.reserve(signal_size);

    // --- Step 4: Precompute constant fitness values from the CSV signal ---
    double model_min, model_max;
    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        model_min = HindmarshRose::MIN;
        model_max = HindmarshRose::MAX;
    }
    else
    {
        throw std::runtime_error("Unsupported model.");
    }

    ConstantSignalFitnessVals stats1 = calc_const_signal_vals(scaled_result.signal, model_min, model_max);

    constexpr size_t POP = GeneticPrivateConfig::POPULATION_SIZE;
    constexpr size_t ELITES = GeneticPrivateConfig::NUM_ELITES;
    constexpr size_t ELITES_SIZE = ELITES * sizeof(Individual);

    // --- Step 5: Initialize population ---
    std::random_device rd;
    std::mt19937 rng(rd());

    std::array<Individual, POP> population = initialize_population<POP>(rng);
    std::array<Individual, POP> new_population;

    // --- Step 6: Evaluate initial population ---
    calc_fitnesses<Integrator, NeuronType, POP>(
        synapsis, model_neur, population, scaled_result,
        stats1, search_phase, model_signal_buffer,
        reset_state_neur, get_v_neur, 0);

    // --- Step 7: Generation loop ---
    std::normal_distribution<double> ndist(0.0, 1.0);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
    std::uniform_real_distribution<double> roulette_dist;

    typename std::array<Individual, POP>::iterator pop_begin = population.begin();
    typename std::array<Individual, POP>::iterator pop_end = population.end();
    typename std::array<Individual, POP>::iterator pop_elites_end = pop_begin + ELITES;

    for (size_t gen = 0; gen < GeneticPrivateConfig::NUM_GENERATIONS; gen++)
    {
        // --- Elitism: use nth_element to get top ELITES ---
        std::nth_element(pop_begin, pop_elites_end, pop_end,
                         fitness_descending);

        // Top ELITES are now at the beginning (not fully sorted)

        // --- Calculate fitness sum for roulette selection ---
        double fitness_sum = 0.0;
        for (const Individual &ind : population)
            fitness_sum += ind.fitness;

        // Assign distribution for this generation
        roulette_dist = std::uniform_real_distribution<double>(0.0, fitness_sum);

        // Keep elites
        std::memcpy(&(new_population[0]), &(population[0]), ELITES_SIZE);

        // Generate offspring
        for (size_t i = ELITES; i < POP; i++)
        {
            const Individual &p1 = roulette_select_one(population, roulette_dist, rng);

            Individual &child = new_population[i];

            // Crossover
            if (prob_dist(rng) < GeneticPrivateConfig::CROSSOVER_PROBABILITY)
            {
                // Select p2 different from p1 using the new function
                const Individual &p2 = roulette_select_second_no_rep(population, p1, roulette_dist, rng);

                // Crossover on params directly
                crossover(p1.params, p2.params, child.params);
            }
            else
            {
                child = p1;
            }

            // Mutation (always attempt, per gene)
            mutate(child.params, rng, ndist, prob_dist);
        }

        std::swap(population, new_population);

        // --- Evaluate only non elites offspring ---
        calc_fitnesses<Integrator, NeuronType, POP>(
            synapsis, model_neur, population, scaled_result,
            stats1, search_phase, model_signal_buffer,
            reset_state_neur, get_v_neur, ELITES);
    }

    // --- Step 8: Return final best individual ---
    return *std::min_element(pop_begin, pop_end, fitness_descending);
}
