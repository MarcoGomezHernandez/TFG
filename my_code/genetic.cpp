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
 * Generate a random ChemicalSynapsisParams within configured ranges
 */
inline ChemicalSynapsisParams random_individual(std::mt19937 &rng)
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

    return p;
}

/*
 * Arithmetic crossover: child = mean of two parents (per parameter), gfast/gslow stay fixed
 */
inline ChemicalSynapsisParams crossover(const ChemicalSynapsisParams &a, const ChemicalSynapsisParams &b)
{
    ChemicalSynapsisParams child;
    child.gfast = GeneticConfig::GFAST_FIXED;
    child.gslow = GeneticConfig::GSLOW_FIXED;
    child.Esyn = (a.Esyn + b.Esyn) / 2.0;
    child.sfast = (a.sfast + b.sfast) / 2.0;
    child.Vfast = (a.Vfast + b.Vfast) / 2.0;
    child.Vslow = (a.Vslow + b.Vslow) / 2.0;
    child.k1 = (a.k1 + b.k1) / 2.0;
    child.k2 = (a.k2 + b.k2) / 2.0;
    child.sslow = (a.sslow + b.sslow) / 2.0;
    return child;
}

/*
 * Mutate an individual with normal perturbation per parameter
 */
inline void mutate(ChemicalSynapsisParams &p, std::mt19937 &rng, std::normal_distribution<double> &ndist)
{
    // gfast and gslow are never mutated
    p.Esyn += ndist(rng) * GeneticConfig::ESYN_MUT_FACTOR;
    p.sfast += ndist(rng) * GeneticConfig::SFAST_MUT_FACTOR;
    p.Vfast += ndist(rng) * GeneticConfig::VFAST_MUT_FACTOR;
    p.Vslow += ndist(rng) * GeneticConfig::VSLOW_MUT_FACTOR;
    p.k1 += ndist(rng) * GeneticConfig::K1_MUT_FACTOR;
    p.k2 += ndist(rng) * GeneticConfig::K2_MUT_FACTOR;
    p.sslow += ndist(rng) * GeneticConfig::SSLOW_MUT_FACTOR;
}

/*
 * Roulette wheel selection: returns index of selected individual
 * fitness_sum is precomputed sum of all fitnesses
 */
inline size_t roulette_select(const std::vector<double> &fitnesses, double fitness_sum, std::mt19937 &rng)
{
    std::uniform_real_distribution<double> dist(0.0, fitness_sum);
    double pick = dist(rng);
    double cumulative = 0.0;
    for (size_t i = 0; i < fitnesses.size(); i++)
    {
        cumulative += fitnesses[i];
        if (cumulative >= pick)
            return i;
    }
    return fitnesses.size() - 1;
}

/*
 * Print a ChemicalSynapsisParams individual
 */
inline void print_individual(const ChemicalSynapsisParams &p, double fitness)
{
    std::cout << "=== Best Individual ===" << std::endl;
    std::cout << "Fitness: " << fitness << std::endl;
    std::cout << "gfast:  " << p.gfast << std::endl;
    std::cout << "Esyn:   " << p.Esyn << std::endl;
    std::cout << "sfast:  " << p.sfast << std::endl;
    std::cout << "Vfast:  " << p.Vfast << std::endl;
    std::cout << "Vslow:  " << p.Vslow << std::endl;
    std::cout << "gslow:  " << p.gslow << std::endl;
    std::cout << "k1:     " << p.k1 << std::endl;
    std::cout << "k2:     " << p.k2 << std::endl;
    std::cout << "sslow:  " << p.sslow << std::endl;
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
void genetic(const std::string &csv_path,
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
        std::cerr << "Error: Signal scaling failed." << std::endl;
        return;
    }

    std::cout << "Signal scaled. dt=" << scaled_result.dt
              << " points_factor=" << scaled_result.points_factor
              << " signal_size=" << scaled_result.signal.size()
              << " interpolated_size=" << scaled_result.interpolated_points.size() << std::endl;

    // --- Step 2: Create neuron and synapsis instances ---
    using SynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    NeuronType model_neur = create_neuron();
    SynapsisType synapsis; // no args, params will be overwritten

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
        std::cerr << "Error: Unsupported model." << std::endl;
        return;
    }

    ConstantSignalFitnessVals stats1 = calc_const_signal_vals(scaled_result.signal, model_min, model_max);

    // --- Step 5: Initialize population ---
    std::random_device rd;
    std::mt19937 rng(rd());

    std::vector<ChemicalSynapsisParams> population(POP);
    std::vector<double> fitnesses(POP, 0.0);

    for (size_t i = 0; i < POP; i++)
    {
        population[i] = random_individual(rng);
    }

    // --- Step 6: Evaluate initial population (all individuals) ---
    // We evaluate in batches of NON_ELITES using calc_fitnesses, but for the first
    // generation we evaluate all POP. We'll do it in a simple loop reusing calc_fitnesses
    // with array size 1 at a time for flexibility, or just batch the whole population.
    // For simplicity, evaluate all POP initially:
    {
        std::array<ChemicalSynapsisParams, POP> params_arr;
        std::array<double, POP> fit_arr;
        for (size_t i = 0; i < POP; i++)
            params_arr[i] = population[i];

        calc_fitnesses<Integrator, NeuronType, POP>(
            synapsis, model_neur, params_arr, scaled_result,
            stats1, search_phase, model_signal_buffer, fit_arr,
            reset_state_neur, get_v_neur);

        for (size_t i = 0; i < POP; i++)
            fitnesses[i] = fit_arr[i];
    }

    // --- Step 7: Generation loop ---
    std::normal_distribution<double> ndist(0.0, 1.0);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

    for (size_t gen = 0; gen < GeneticConfig::NUM_GENERATIONS; gen++)
    {
        // --- Elitism: find top ELITES individuals ---
        // Get indices sorted by fitness descending
        std::vector<size_t> sorted_indices(POP);
        std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
        std::sort(sorted_indices.begin(), sorted_indices.end(),
                  [&fitnesses](size_t a, size_t b)
                  { return fitnesses[a] > fitnesses[b]; });

        // Save elites (params and fitnesses)
        std::array<ChemicalSynapsisParams, ELITES> elite_params;
        std::array<double, ELITES> elite_fitnesses;
        for (size_t e = 0; e < ELITES; e++)
        {
            elite_params[e] = population[sorted_indices[e]];
            elite_fitnesses[e] = fitnesses[sorted_indices[e]];
        }

        // --- Selection, crossover, mutation to produce NON_ELITES offspring ---
        double fitness_sum = 0.0;
        for (size_t i = 0; i < POP; i++)
            fitness_sum += fitnesses[i];

        std::vector<ChemicalSynapsisParams> new_population(POP);

        // Place elites first (they keep their params and fitness, no recalculation)
        for (size_t e = 0; e < ELITES; e++)
        {
            new_population[e] = elite_params[e];
        }

        // Generate NON_ELITES offspring
        for (size_t i = ELITES; i < POP; i++)
        {
            // Roulette selection of two parents (with repetition, from full population including elites)
            size_t p1 = roulette_select(fitnesses, fitness_sum, rng);
            size_t p2 = roulette_select(fitnesses, fitness_sum, rng);

            ChemicalSynapsisParams child;

            // Crossover
            if (prob_dist(rng) < GeneticConfig::CROSSOVER_PROBABILITY)
            {
                child = crossover(population[p1], population[p2]);
            }
            else
            {
                // No crossover: copy one parent
                child = population[p1];
            }

            // Mutation
            if (prob_dist(rng) < GeneticConfig::MUTATION_PROBABILITY)
            {
                mutate(child, rng, ndist);
            }

            // Ensure gfast and gslow stay fixed
            child.gfast = GeneticConfig::GFAST_FIXED;
            child.gslow = GeneticConfig::GSLOW_FIXED;

            new_population[i] = child;
        }

        population = new_population;

        // --- Evaluate only NON_ELITES offspring ---
        {
            std::array<ChemicalSynapsisParams, NON_ELITES> ne_params;
            std::array<double, NON_ELITES> ne_fitnesses;
            for (size_t i = 0; i < NON_ELITES; i++)
                ne_params[i] = population[ELITES + i];

            calc_fitnesses<Integrator, NeuronType, NON_ELITES>(
                synapsis, model_neur, ne_params, scaled_result,
                stats1, search_phase, model_signal_buffer, ne_fitnesses,
                reset_state_neur, get_v_neur);

            for (size_t i = 0; i < NON_ELITES; i++)
                fitnesses[ELITES + i] = ne_fitnesses[i];
        }

        // Restore elite fitnesses (not recalculated)
        for (size_t e = 0; e < ELITES; e++)
        {
            fitnesses[e] = elite_fitnesses[e];
        }

        // Print best fitness of this generation
        double best_fit = *std::max_element(fitnesses.begin(), fitnesses.end());
        std::cout << "Generation " << gen + 1 << "/" << GeneticConfig::NUM_GENERATIONS
                  << " | Best fitness: " << best_fit << std::endl;
    }

    // --- Step 8: Print final best individual ---
    size_t best_idx = std::distance(fitnesses.begin(),
                                    std::max_element(fitnesses.begin(), fitnesses.end()));
    print_individual(population[best_idx], fitnesses[best_idx]);
}
