#include <array>
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <utility>
#include <stdexcept>
#include <cstring>
#include "scaling.hpp"
#include <ChemicalSynapsis.h>

/*
 * Genetic algorithm configuration constants (private)
 */
namespace GeneticConfig
{
    // Population and generations
    static constexpr size_t POPULATION_SIZE = 40;
    static constexpr size_t NUM_GENERATIONS = 80;
    static constexpr size_t NUM_ELITES = 3;

    // Observation time is use_time / this divisor
    static constexpr double OBSERVATION_TIME_DIVISOR = 3.0;

    // Crossover and mutation probabilities
    static constexpr double CROSSOVER_PROBABILITY = 0.9;
    static constexpr double MUTATION_PROBABILITY = 0.1;

    // Mutation scale factor (η): percentage of parameter range
    static constexpr double ETA = 0.2;

    // Random initialization ranges [min, max] for each mutable parameter
    // ESyn tiene rangos distintos según estemos buscando en fase (search_phase = true)
    // o en antifase (search_phase = false). PDF usa -1.92 para red asimétrica.
    static constexpr double ESYN_MIN_PHASE = 0.0;
    static constexpr double ESYN_MAX_PHASE = 5.0;
    static constexpr double ESYN_MIN_ANTIPHASE = -5.0;
    static constexpr double ESYN_MAX_ANTIPHASE = -0.5;

    static constexpr double SFAST_MIN = 0.01;
    static constexpr double SFAST_MAX = 10.0;

    static constexpr double SSLOW_MIN = 0.01;
    static constexpr double SSLOW_MAX = 10.0;

    static constexpr double VFAST_MIN = -3.0;
    static constexpr double VFAST_MAX = 3.0;

    static constexpr double VSLOW_MIN = -3.0;
    static constexpr double VSLOW_MAX = 3.0;

    static constexpr double K1_MIN = 0.001;
    static constexpr double K1_MAX = 10.0;

    static constexpr double K2_MIN = 0.0001;
    static constexpr double K2_MAX = 1.0;

    static constexpr double GFAST_MIN = 1e-3;
    static constexpr double GFAST_MAX = 5.0;

    static constexpr double GSLOW_MIN = 1e-3;
    static constexpr double GSLOW_MAX = 5.0;

    // Mutation perturbation factors: σ_p = η × (p_max - p_min)
    static constexpr double ESYN_MUT_FACTOR_PHASE = ETA * (ESYN_MAX_PHASE - ESYN_MIN_PHASE);
    static constexpr double ESYN_MUT_FACTOR_ANTIPHASE = ETA * (ESYN_MAX_ANTIPHASE - ESYN_MIN_ANTIPHASE);
    static constexpr double SFAST_MUT_FACTOR = ETA * (SFAST_MAX - SFAST_MIN);
    static constexpr double VFAST_MUT_FACTOR = ETA * (VFAST_MAX - VFAST_MIN);
    static constexpr double VSLOW_MUT_FACTOR = ETA * (VSLOW_MAX - VSLOW_MIN);
    static constexpr double K1_MUT_FACTOR = ETA * (K1_MAX - K1_MIN);
    static constexpr double K2_MUT_FACTOR = ETA * (K2_MAX - K2_MIN);
    static constexpr double SSLOW_MUT_FACTOR = ETA * (SSLOW_MAX - SSLOW_MIN);
    static constexpr double GFAST_MUT_FACTOR = ETA * (GFAST_MAX - GFAST_MIN);
    static constexpr double GSLOW_MUT_FACTOR = ETA * (GSLOW_MAX - GSLOW_MIN);
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
inline std::array<Individual, POP_SIZE> initialize_population(std::mt19937 &rng, bool search_phase, SynComponent syn_component)
{
    double esyn_min, esyn_max;
    if (search_phase)
    {
        esyn_min = GeneticConfig::ESYN_MIN_PHASE;
        esyn_max = GeneticConfig::ESYN_MAX_PHASE;
    }
    else
    {
        esyn_min = GeneticConfig::ESYN_MIN_ANTIPHASE;
        esyn_max = GeneticConfig::ESYN_MAX_ANTIPHASE;
    }
    // Always needed
    std::uniform_real_distribution<double> dist_esyn(esyn_min, esyn_max);

    // Declare all distributions (required), but only initialize those we'll use.
    std::uniform_real_distribution<double> dist_sfast;
    std::uniform_real_distribution<double> dist_vfast;
    std::uniform_real_distribution<double> dist_vslow;
    std::uniform_real_distribution<double> dist_k1;
    std::uniform_real_distribution<double> dist_k2;
    std::uniform_real_distribution<double> dist_sslow;
    std::uniform_real_distribution<double> dist_gfast;
    std::uniform_real_distribution<double> dist_gslow;

    // Initialize only the distributions required by the requested syn_component
    if (syn_component != SynComponent::ISLOW) // using IFAST or BOTH
    {
        dist_sfast = std::uniform_real_distribution<double>(GeneticConfig::SFAST_MIN, GeneticConfig::SFAST_MAX);
        dist_vfast = std::uniform_real_distribution<double>(GeneticConfig::VFAST_MIN, GeneticConfig::VFAST_MAX);
        dist_gfast = std::uniform_real_distribution<double>(GeneticConfig::GFAST_MIN, GeneticConfig::GFAST_MAX);
    }

    if (syn_component != SynComponent::IFAST) // using ISLOW or BOTH
    {
        dist_vslow = std::uniform_real_distribution<double>(GeneticConfig::VSLOW_MIN, GeneticConfig::VSLOW_MAX);
        dist_k1 = std::uniform_real_distribution<double>(GeneticConfig::K1_MIN, GeneticConfig::K1_MAX);
        dist_k2 = std::uniform_real_distribution<double>(GeneticConfig::K2_MIN, GeneticConfig::K2_MAX);
        dist_sslow = std::uniform_real_distribution<double>(GeneticConfig::SSLOW_MIN, GeneticConfig::SSLOW_MAX);
        dist_gslow = std::uniform_real_distribution<double>(GeneticConfig::GSLOW_MIN, GeneticConfig::GSLOW_MAX);
    }

    const bool use_ifast = (syn_component != SynComponent::ISLOW);
    const bool use_islow = (syn_component != SynComponent::IFAST);

    std::array<Individual, POP_SIZE> population;

    for (Individual &ind : population)
    {
        ChemicalSynapsisParams &p = ind.params;
        p.Esyn = dist_esyn(rng);
        if (use_ifast)
        {
            p.gfast = dist_gfast(rng);
            p.sfast = dist_sfast(rng);
            p.Vfast = dist_vfast(rng);
        }
        else
        {
            p.gfast = 0.0;
            p.sfast = 0.0;
            p.Vfast = 0.0;
        }

        if (use_islow)
        {
            p.gslow = dist_gslow(rng);
            p.Vslow = dist_vslow(rng);
            p.k1 = dist_k1(rng);
            p.k2 = dist_k2(rng);
            p.sslow = dist_sslow(rng);
        }
        else
        {
            p.gslow = 0.0;
            p.Vslow = 0.0;
            p.k1 = 0.0;
            p.k2 = 0.0;
            p.sslow = 0.0;
        }
    }

    return population;
}

/*
 * Arithmetic crossover: operate directly on ChemicalSynapsisParams (child = mean of two parents)
 */
inline void crossover(const ChemicalSynapsisParams &a, const ChemicalSynapsisParams &b, ChemicalSynapsisParams &result)
{
    result.gfast = (a.gfast + b.gfast) * 0.5;
    result.gslow = (a.gslow + b.gslow) * 0.5;
    result.Esyn = (a.Esyn + b.Esyn) * 0.5;
    result.sfast = (a.sfast + b.sfast) * 0.5;
    result.Vfast = (a.Vfast + b.Vfast) * 0.5;
    result.Vslow = (a.Vslow + b.Vslow) * 0.5;
    result.k1 = (a.k1 + b.k1) * 0.5;
    result.k2 = (a.k2 + b.k2) * 0.5;
    result.sslow = (a.sslow + b.sslow) * 0.5;
}

/*
 * Mutate: recibe directamente los parámetros de la sinapsis (params) y aplica
 * perturbaciones normales por cada campo con su probabilidad independiente.
 */
inline void mutate(ChemicalSynapsisParams &p, std::mt19937 &rng, std::normal_distribution<double> &ndist, std::uniform_real_distribution<double> &prob_dist, double esyn_mut_factor, SynComponent syn_component)
{
    constexpr double MUTATION_PROBABILITY = GeneticConfig::MUTATION_PROBABILITY;
    // Precomputed local aliases for GeneticConfig minima used multiple times
    constexpr double GFAST_MIN = GeneticConfig::GFAST_MIN;
    constexpr double SFAST_MIN = GeneticConfig::SFAST_MIN;
    constexpr double GSLOW_MIN = GeneticConfig::GSLOW_MIN;
    constexpr double K1_MIN = GeneticConfig::K1_MIN;
    constexpr double K2_MIN = GeneticConfig::K2_MIN;
    constexpr double SSLOW_MIN = GeneticConfig::SSLOW_MIN;

    // Esyn is shared by both components
    if (prob_dist(rng) < MUTATION_PROBABILITY)
        p.Esyn += ndist(rng) * esyn_mut_factor;
    // ifast params
    if (syn_component != SynComponent::ISLOW)
    {
        double &gfast = p.gfast;
        double &sfast = p.sfast;

        if (prob_dist(rng) < MUTATION_PROBABILITY)
            gfast += ndist(rng) * GeneticConfig::GFAST_MUT_FACTOR;
        if (prob_dist(rng) < MUTATION_PROBABILITY)
            sfast += ndist(rng) * GeneticConfig::SFAST_MUT_FACTOR;
        if (prob_dist(rng) < MUTATION_PROBABILITY)
            p.Vfast += ndist(rng) * GeneticConfig::VFAST_MUT_FACTOR;

        if (gfast < GFAST_MIN)
            gfast = GFAST_MIN;
        if (sfast < SFAST_MIN)
            sfast = SFAST_MIN;
    }
    // islow params
    if (syn_component != SynComponent::IFAST)
    {
        double &gslow = p.gslow;
        double &k1 = p.k1;
        double &k2 = p.k2;
        double &sslow = p.sslow;

        if (prob_dist(rng) < MUTATION_PROBABILITY)
            gslow += ndist(rng) * GeneticConfig::GSLOW_MUT_FACTOR;
        if (prob_dist(rng) < MUTATION_PROBABILITY)
            p.Vslow += ndist(rng) * GeneticConfig::VSLOW_MUT_FACTOR;
        if (prob_dist(rng) < MUTATION_PROBABILITY)
            k1 += ndist(rng) * GeneticConfig::K1_MUT_FACTOR;
        if (prob_dist(rng) < MUTATION_PROBABILITY)
            k2 += ndist(rng) * GeneticConfig::K2_MUT_FACTOR;
        if (prob_dist(rng) < MUTATION_PROBABILITY)
            sslow += ndist(rng) * GeneticConfig::SSLOW_MUT_FACTOR;

        if (gslow < GSLOW_MIN)
            gslow = GSLOW_MIN;
        if (k1 < K1_MIN)
            k1 = K1_MIN;
        if (k2 < K2_MIN)
            k2 = K2_MIN;
        if (sslow < SSLOW_MIN)
            sslow = SSLOW_MIN;
    }
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
    const double pick = roulette_dist(rng);
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
 * Uses a distribution with fitness sum excluding the ignored individual
 */
template <size_t N>
inline const Individual &roulette_select_second_no_rep(const std::array<Individual, N> &population,
                                                       const Individual &ignore_ind,
                                                       std::uniform_real_distribution<double> &roulette_dist,
                                                       std::mt19937 &rng)
{
    const double pick = roulette_dist(rng);
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
                   SynComponent syn_component,
                   CreateFuncType create_neur,
                   ResetStateFuncType reset_state_neur,
                   GetVFuncType get_v_neur,
                   typename NeuronType::variable neur_v_var,
                   int syn_model_step_factor,
                   bool verbose)
{
    // Calculate observation time
    const double observation_time = use_time / GeneticConfig::OBSERVATION_TIME_DIVISOR;

    // --- Step 1: Scale the signal ---
    ScaledSignalResult scaled_result = scale_signal(
        csv_path, column_index, csv_step, start_time, use_time,
        observation_time, integrator, model, check_drift);

    if (!scaled_result.success)
    {
        throw std::runtime_error("Signal scaling failed.");
    }

    // --- Step 2: Create neuron and synapsis instances ---
    NeuronType model_neur = create_neur(false);
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;
    typename ChemicalSynapsisType::ConstructorArgs syn_args{};
    ChemicalSynapsisType synapsis(create_neur(true), neur_v_var, model_neur, neur_v_var, syn_args, syn_model_step_factor);

    // --- Step 3: Allocate buffers ---
    const size_t signal_size = scaled_result.signal.size();
    std::vector<double> model_signal_buffer(signal_size);
    std::vector<double> synapsis_signal_buffer(signal_size);

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

    const size_t avg_smooth_points_living = std::max(FitnessConstants::MIN_AVG_SMOOTH_POINTS, static_cast<size_t>(scaled_result.pts_burst_real * FitnessConfig::AVG_SMOOTH_POINTS_BURST_FRACTION));
    const size_t avg_smooth_points_model = avg_smoothing_points_living * scaled_result.points_factor;

    const ConstantSignalFitnessVals living_const_signal_fitness_vals = calc_const_signal_fitness_vals(scaled_result.signal, model_min, model_max, search_phase,
                                                                                                      avg_smooth_points_living);

    constexpr size_t POP = GeneticConfig::POPULATION_SIZE;
    constexpr size_t ELITES = GeneticConfig::NUM_ELITES;
    constexpr size_t ELITES_SIZE = ELITES * sizeof(Individual);

    // --- Step 5: Initialize population ---
    std::random_device rd;
    std::mt19937 rng(rd());

    double esyn_mut_factor;
    if (search_phase)
    {
        esyn_mut_factor = GeneticConfig::ESYN_MUT_FACTOR_PHASE;
    }
    else
    {
        esyn_mut_factor = GeneticConfig::ESYN_MUT_FACTOR_ANTIPHASE;
    }

    std::array<Individual, POP> buffer1 = initialize_population<POP>(rng, search_phase, syn_component);
    std::array<Individual, POP> buffer2;

    std::array<Individual, POP> *population_ptr = &buffer1;
    std::array<Individual, POP> *new_population_ptr = &buffer2;

    // --- Step 6: Evaluate initial population ---
    calc_fitnesses<Integrator, NeuronType, POP>(
        synapsis, model_neur, *population_ptr, scaled_result,
        living_const_signal_fitness_vals, search_phase, model_signal_buffer, synapsis_signal_buffer,
        reset_state_neur, get_v_neur, 0, avg_smooth_points_model);

    // --- Step 7: Generation loop ---
    std::normal_distribution<double> ndist(0.0, 1.0);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
    std::uniform_real_distribution<double> roulette_dist;
    std::uniform_real_distribution<double> roulette_dist_no_p1;

    typename std::array<Individual, POP>::iterator pop_begin;

    for (size_t gen = 0; gen < GeneticConfig::NUM_GENERATIONS; gen++)
    {
        if (verbose)
        {
            std::cout << "\rGeneration " << (gen + 1) << "/" << GeneticConfig::NUM_GENERATIONS << std::flush;
        }

        std::array<Individual, POP> &population = *population_ptr;
        std::array<Individual, POP> &new_population = *new_population_ptr;

        // --- Elitism: use nth_element to get top ELITES ---
        pop_begin = population.begin();
        std::nth_element(pop_begin, pop_begin + ELITES, population.end(),
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
            if (prob_dist(rng) < GeneticConfig::CROSSOVER_PROBABILITY)
            {
                // Calculate fitness sum excluding p1 and assign new distribution
                // Select p2 different from p1
                roulette_dist_no_p1 = std::uniform_real_distribution<double>(0.0, fitness_sum - p1.fitness);
                const Individual &p2 = roulette_select_second_no_rep(population, p1, roulette_dist_no_p1, rng);

                // Crossover on params directly
                crossover(p1.params, p2.params, child.params);
            }
            else
            {
                child = p1;
            }

            // Mutation (always attempt, per gene)
            mutate(child.params, rng, ndist, prob_dist, esyn_mut_factor, syn_component);
        }

        std::swap(population_ptr, new_population_ptr);

        // --- Evaluate only non elites offspring ---
        calc_fitnesses<Integrator, NeuronType, POP>(
            synapsis, model_neur, *population_ptr, scaled_result,
            living_const_signal_fitness_vals, search_phase, model_signal_buffer, synapsis_signal_buffer,
            reset_state_neur, get_v_neur, ELITES, avg_smooth_points_model);
    }

    // --- Step 8: Return final best individual ---
    std::array<Individual, POP> &population = *population_ptr;
    const Individual &best = *std::min_element(population.begin(), population.end(), fitness_descending);

    if (verbose)
    {
        std::cout << std::endl;
        const ChemicalSynapsisParams &params = best.params;
        std::cout << "Fitness: " << best.fitness << std::endl;
        std::cout << "gfast_values = [" << params.gfast << "]" << std::endl;
        std::cout << "gslow_values = [" << params.gslow << "]" << std::endl;
        std::cout << "Esyn_values = [" << params.Esyn << "]" << std::endl;
        std::cout << "Vfast_values = [" << params.Vfast << "]" << std::endl;
        std::cout << "Vslow_values = [" << params.Vslow << "]" << std::endl;
        std::cout << "sfast_values = [" << params.sfast << "]" << std::endl;
        std::cout << "sslow_values = [" << params.sslow << "]" << std::endl;
        std::cout << "k1_values = [" << params.k1 << "]" << std::endl;
        std::cout << "k2_values = [" << params.k2 << "]" << std::endl;
    }

    return best;
}
