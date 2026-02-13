#include "genetic.hpp"

/*
 * Genetic algorithm configuration constants
 */
namespace GeneticConfig
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

    // Fixed values for gfast and gslow
    static constexpr double GFAST_FIXED = 0.1;
    static constexpr double GSLOW_FIXED = 0.01;

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
 * Generate a random Individual within configured ranges
 */
inline Individual random_individual(std::mt19937 &rng)
{
    ChemicalSynapsisParams p;
    p.gfast = GeneticConfig::GFAST_FIXED;
    p.gslow = GeneticConfig::GSLOW_FIXED;

    std::uniform_real_distribution<double> dist_esyn(GeneticConfig::ESYN_MIN, GeneticConfig::ESYN_MAX);
    std::uniform_real_distribution<double> dist_sfast(GeneticConfig::SFAST_MIN, GeneticConfig::SFAST_MAX);
    std::uniform_real_distribution<double> dist_vfast(GeneticConfig::VFAST_MIN, GeneticConfig::VFAST_MAX);
    std::uniform_real_distribution<double> dist_vslow(GeneticConfig::VSLOW_MIN, GeneticConfig::VSLOW_MAX);
    std::uniform_real_distribution<double> dist_k1(GeneticConfig::K1_MIN, GeneticConfig::K1_MAX);
    std::uniform_real_distribution<double> dist_k2(GeneticConfig::K2_MIN, GeneticConfig::K2_MAX);
    std::uniform_real_distribution<double> dist_sslow(GeneticConfig::SSLOW_MIN, GeneticConfig::SSLOW_MAX);

    p.Esyn = dist_esyn(rng);
    p.sfast = dist_sfast(rng);
    p.Vfast = dist_vfast(rng);
    p.Vslow = dist_vslow(rng);
    p.k1 = dist_k1(rng);
    p.k2 = dist_k2(rng);
    p.sslow = dist_sslow(rng);

    return {p, 0.0};
}

/*
 * Arithmetic crossover: child = mean of two parents (per parameter), gfast/gslow stay fixed
 */
inline Individual crossover(const Individual &a, const Individual &b)
{
    ChemicalSynapsisParams child;
    child.gfast = GeneticConfig::GFAST_FIXED;
    child.gslow = GeneticConfig::GSLOW_FIXED;
    child.Esyn = (a.params.Esyn + b.params.Esyn) / 2.0;
    child.sfast = (a.params.sfast + b.params.sfast) / 2.0;
    child.Vfast = (a.params.Vfast + b.params.Vfast) / 2.0;
    child.Vslow = (a.params.Vslow + b.params.Vslow) / 2.0;
    child.k1 = (a.params.k1 + b.params.k1) / 2.0;
    child.k2 = (a.params.k2 + b.params.k2) / 2.0;
    child.sslow = (a.params.sslow + b.params.sslow) / 2.0;
    return {child, 0.0};
}

/*
 * Mutate an individual with normal perturbation per parameter
 */
inline void mutate(Individual &ind, std::mt19937 &rng, std::normal_distribution<double> &ndist)
{
    // gfast and gslow are never mutated
    ind.params.Esyn += ndist(rng) * GeneticConfig::ESYN_MUT_FACTOR;
    ind.params.sfast += ndist(rng) * GeneticConfig::SFAST_MUT_FACTOR;
    ind.params.Vfast += ndist(rng) * GeneticConfig::VFAST_MUT_FACTOR;
    ind.params.Vslow += ndist(rng) * GeneticConfig::VSLOW_MUT_FACTOR;
    ind.params.k1 += ndist(rng) * GeneticConfig::K1_MUT_FACTOR;
    ind.params.k2 += ndist(rng) * GeneticConfig::K2_MUT_FACTOR;
    ind.params.sslow += ndist(rng) * GeneticConfig::SSLOW_MUT_FACTOR;
}

/*
 * Roulette wheel selection: returns index of selected individual
 * fitness_sum is precomputed sum of all fitnesses
 */
template <size_t N>
inline size_t roulette_select(const std::array<Individual, N> &population, double fitness_sum, std::mt19937 &rng)
{
    std::uniform_real_distribution<double> dist(0.0, fitness_sum);
    double pick = dist(rng);
    double cumulative = 0.0;
    for (size_t i = 0; i < N; i++)
    {
        cumulative += population[i].fitness;
        if (cumulative >= pick)
            return i;
    }
    return N - 1;
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
                   NumericIntegrator integrator_enum,
                   NeuronModel model,
                   bool search_phase,
                   bool check_drift,
                   CreateFuncType create_neuron,
                   ResetStateFuncType reset_state_neur,
                   GetVFuncType get_v_neur)
{
    constexpr size_t POP = GeneticConfig::POPULATION_SIZE;
    constexpr size_t ELITES = GeneticConfig::NUM_ELITES;
    constexpr size_t NON_ELITES = POP - ELITES;

    // Calculate observation time
    const double observation_time = use_time / GeneticConfig::OBSERVATION_TIME_DIVISOR;

    // --- Step 1: Scale the signal ---
    ScaledSignalResult scaled_result = scale_signal(
        csv_path, column_index, csv_step, start_time, use_time,
        observation_time, integrator_enum, model, check_drift);

    if (!scaled_result.success)
    {
        throw std::runtime_error("Signal scaling failed.");
    }

    // --- Step 2: Create neuron and synapsis instances ---
    using SynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    NeuronType model_neur = create_neuron();
    SynapsisType synapsis;

    // --- Step 3: Allocate buffers ---
    const size_t signal_size = scaled_result.signal.size();
    std::vector<double> model_signal_buffer(signal_size);

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

    // --- Step 5: Initialize population ---
    std::random_device rd;
    std::mt19937 rng(rd());

    std::array<Individual, POP> population;

    for (size_t i = 0; i < POP; i++)
    {
        population[i] = random_individual(rng);
    }

    // --- Step 6: Evaluate initial population ---
    calc_fitnesses<Integrator, NeuronType, POP>(
        synapsis, model_neur, population, scaled_result,
        stats1, search_phase, model_signal_buffer,
        reset_state_neur, get_v_neur);

    // --- Step 7: Generation loop ---
    std::normal_distribution<double> ndist(0.0, 1.0);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

    for (size_t gen = 0; gen < GeneticConfig::NUM_GENERATIONS; gen++)
    {
        // --- Elitism: sort population by fitness descending ---
        std::sort(population.begin(), population.end(),
                  [](const Individual &a, const Individual &b)
                  { return a.fitness > b.fitness; });

        // Top ELITES are already at the beginning after sorting

        // --- Calculate fitness sum for roulette selection ---
        double fitness_sum = 0.0;
        for (const auto &ind : population)
            fitness_sum += ind.fitness;

        // --- Generate NON_ELITES offspring ---
        std::array<Individual, POP> new_population;

        // Keep elites
        for (size_t e = 0; e < ELITES; e++)
        {
            new_population[e] = population[e];
        }

        // Generate offspring
        for (size_t i = ELITES; i < POP; i++)
        {
            size_t p1 = roulette_select(population, fitness_sum, rng);
            size_t p2 = roulette_select(population, fitness_sum, rng);

            Individual child;

            // Crossover
            if (prob_dist(rng) < GeneticConfig::CROSSOVER_PROBABILITY)
            {
                child = crossover(population[p1], population[p2]);
            }
            else
            {
                child = population[p1];
            }

            // Mutation
            if (prob_dist(rng) < GeneticConfig::MUTATION_PROBABILITY)
            {
                mutate(child, rng, ndist);
            }

            // Ensure gfast and gslow stay fixed
            child.params.gfast = GeneticConfig::GFAST_FIXED;
            child.params.gslow = GeneticConfig::GSLOW_FIXED;

            new_population[i] = child;
        }

        population = new_population;

        // --- Evaluate only NON_ELITES offspring ---
        std::array<Individual, NON_ELITES> ne_arr;
        for (size_t i = 0; i < NON_ELITES; i++)
            ne_arr[i] = population[ELITES + i];

        calc_fitnesses<Integrator, NeuronType, NON_ELITES>(
            synapsis, model_neur, ne_arr, scaled_result,
            stats1, search_phase, model_signal_buffer,
            reset_state_neur, get_v_neur);

        for (size_t i = 0; i < NON_ELITES; i++)
            population[ELITES + i] = ne_arr[i];
    }

    // --- Step 8: Return final best individual ---
    auto best_it = std::max_element(population.begin(), population.end(),
                                    [](const Individual &a, const Individual &b)
                                    { return a.fitness < b.fitness; });
    return *best_it;
}
