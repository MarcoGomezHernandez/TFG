#include <array>
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

namespace GeneticConfig
{
    static constexpr size_t POPULATION_SIZE = 40;
    static constexpr size_t NUM_GENERATIONS = 80;
    static constexpr size_t NUM_ELITES = 3;

    static constexpr double OBSERVATION_TIME_DIVISOR = 3.0;

    static constexpr double CROSSOVER_PROBABILITY = 0.9;
    static constexpr double MUTATION_PROBABILITY = 0.1;

    static constexpr double ETA = 0.2;

    static constexpr double ESYN_MIN_PHASE = -5.0;
    static constexpr double ESYN_MAX_PHASE = 5.0;
    static constexpr double ESYN_MIN_ANTIPHASE = -5.0;
    static constexpr double ESYN_MAX_ANTIPHASE = 5.0;

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

static bool fitness_descending(const Individual &a, const Individual &b)
{
    return a.fitness > b.fitness;
}

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
    std::uniform_real_distribution<double> dist_esyn(esyn_min, esyn_max);

    std::uniform_real_distribution<double> dist_sfast;
    std::uniform_real_distribution<double> dist_vfast;
    std::uniform_real_distribution<double> dist_vslow;
    std::uniform_real_distribution<double> dist_k1;
    std::uniform_real_distribution<double> dist_k2;
    std::uniform_real_distribution<double> dist_sslow;
    std::uniform_real_distribution<double> dist_gfast;
    std::uniform_real_distribution<double> dist_gslow;

    if (syn_component != SynComponent::ISLOW)
    {
        dist_sfast = std::uniform_real_distribution<double>(GeneticConfig::SFAST_MIN, GeneticConfig::SFAST_MAX);
        dist_vfast = std::uniform_real_distribution<double>(GeneticConfig::VFAST_MIN, GeneticConfig::VFAST_MAX);
        dist_gfast = std::uniform_real_distribution<double>(GeneticConfig::GFAST_MIN, GeneticConfig::GFAST_MAX);
    }

    if (syn_component != SynComponent::IFAST)
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

inline void mutate(ChemicalSynapsisParams &p, std::mt19937 &rng, std::normal_distribution<double> &ndist, std::uniform_real_distribution<double> &prob_dist, double esyn_mut_factor, bool use_ifast, bool use_islow)
{
    constexpr double MUTATION_PROBABILITY = GeneticConfig::MUTATION_PROBABILITY;
    constexpr double GFAST_MIN = GeneticConfig::GFAST_MIN;
    constexpr double SFAST_MIN = GeneticConfig::SFAST_MIN;
    constexpr double GSLOW_MIN = GeneticConfig::GSLOW_MIN;
    constexpr double K1_MIN = GeneticConfig::K1_MIN;
    constexpr double K2_MIN = GeneticConfig::K2_MIN;
    constexpr double SSLOW_MIN = GeneticConfig::SSLOW_MIN;

    if (prob_dist(rng) < MUTATION_PROBABILITY)
        p.Esyn += ndist(rng) * esyn_mut_factor;
    if (use_ifast)
    {
        double &gfast = p.gfast;
        double &sfast = p.sfast;

        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double delta = ndist(rng) * GeneticConfig::GFAST_MUT_FACTOR;
            if (gfast == GFAST_MIN && delta < 0)
                gfast -= delta;
            else
            {
                gfast += delta;
                if (gfast < GFAST_MIN)
                    gfast = GFAST_MIN;
            }
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double delta = ndist(rng) * GeneticConfig::SFAST_MUT_FACTOR;
            if (sfast == SFAST_MIN && delta < 0)
                sfast -= delta;
            else
            {
                sfast += delta;
                if (sfast < SFAST_MIN)
                    sfast = SFAST_MIN;
            }
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
            p.Vfast += ndist(rng) * GeneticConfig::VFAST_MUT_FACTOR;
    }
    if (use_islow)
    {
        double &gslow = p.gslow;
        double &k1 = p.k1;
        double &k2 = p.k2;
        double &sslow = p.sslow;

        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double delta = ndist(rng) * GeneticConfig::GSLOW_MUT_FACTOR;
            if (gslow == GSLOW_MIN && delta < 0)
                gslow -= delta;
            else
            {
                gslow += delta;
                if (gslow < GSLOW_MIN)
                    gslow = GSLOW_MIN;
            }
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
            p.Vslow += ndist(rng) * GeneticConfig::VSLOW_MUT_FACTOR;
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double delta = ndist(rng) * GeneticConfig::K1_MUT_FACTOR;
            if (k1 == K1_MIN && delta < 0)
                k1 -= delta;
            else
            {
                k1 += delta;
                if (k1 < K1_MIN)
                    k1 = K1_MIN;
            }
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double delta = ndist(rng) * GeneticConfig::K2_MUT_FACTOR;
            if (k2 == K2_MIN && delta < 0)
                k2 -= delta;
            else
            {
                k2 += delta;
                if (k2 < K2_MIN)
                    k2 = K2_MIN;
            }
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double delta = ndist(rng) * GeneticConfig::SSLOW_MUT_FACTOR;
            if (sslow == SSLOW_MIN && delta < 0)
                sslow -= delta;
            else
            {
                sslow += delta;
                if (sslow < SSLOW_MIN)
                    sslow = SSLOW_MIN;
            }
        }
    }
}

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
    return population.back();
}

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

    for (const Individual &ind : population)
    {
        if (&ind != &ignore_ind)
            return ind;
    }
    return population.back();
}

template <typename Integrator, typename NeuronType,
          CreateFunc<NeuronType> CreateFuncType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
Individual genetic(const std::string &csv_path,
                   size_t column_i,
                   double csv_step,
                   double start_time,
                   double use_time,
                   double stabilization_time,
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
    const double observation_time = use_time / GeneticConfig::OBSERVATION_TIME_DIVISOR;

    ScaledSigResult scaled_result = scale_sig(
        csv_path, column_i, csv_step, start_time, use_time + stabilization_time,
        observation_time, integrator, model, check_drift);

    if (!scaled_result.success)
    {
        throw std::runtime_error("Signal scaling failed.");
    }

    NeuronType model_neur = create_neur(false);
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;
    typename ChemicalSynapsisType::ConstructorArgs syn_args{};
    ChemicalSynapsisType synapsis(create_neur(true), neur_v_var, model_neur, neur_v_var, syn_args, syn_model_step_factor);

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

    const size_t avg_smooth_points = std::max(FitnessConstants::MIN_AVG_SMOOTH_POINTS, static_cast<size_t>(scaled_result.pts_burst_real / FitnessConfig::AVG_SMOOTH_POINTS_BURST_DIVISOR));

    const size_t total_vpre_sig_size = scaled_result.vpre_sig.size();
    const size_t stabilization_points = std::max(static_cast<size_t>(stabilization_time / csv_step), avg_smooth_points);
    const size_t use_vpre_sig_size = total_vpre_sig_size - stabilization_points;

    SigBuffers buffers;
    buffers.vpost_sig.resize(avg_smooth_points + use_vpre_sig_size);
    buffers.i_sig.resize(use_vpre_sig_size);
    const bool use_ifast = (syn_component != SynComponent::ISLOW);
    const bool use_islow = (syn_component != SynComponent::IFAST);
    if (use_ifast)
        buffers.ifast_sig.resize(use_vpre_sig_size);
    if (use_islow)
        buffers.islow_sig.resize(use_vpre_sig_size);
    buffers.kfr_padded.resize(use_vpre_sig_size + 2 * FitnessConfig::FILTER_PAD_LEN);

    const ConstantSigFitnessVals living_const_sig_fitness_vals = calc_const_sig_fitness_vals(
        scaled_result.vpre_sig, model_min, model_max, search_phase,
        avg_smooth_points, stabilization_points, scaled_result.pts_burst_real, buffers, use_ifast, use_islow);

    constexpr size_t POP = GeneticConfig::POPULATION_SIZE;
    constexpr size_t ELITES = GeneticConfig::NUM_ELITES;
    constexpr size_t ELITES_SIZE = ELITES * sizeof(Individual);

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

    calc_fitnesses<Integrator, NeuronType, POP>(
        synapsis, model_neur, *population_ptr, scaled_result,
        living_const_sig_fitness_vals, search_phase, buffers,
        reset_state_neur, get_v_neur, 0, stabilization_points, avg_smooth_points,
        use_ifast, use_islow);

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

        pop_begin = population.begin();
        std::nth_element(pop_begin, pop_begin + ELITES, population.end(),
                         fitness_descending);

        double fitness_sum = 0.0;
        for (const Individual &ind : population)
            fitness_sum += ind.fitness;

        roulette_dist = std::uniform_real_distribution<double>(0.0, fitness_sum);

        std::memcpy(&(new_population[0]), &(population[0]), ELITES_SIZE);

        for (size_t i = ELITES; i < POP; i++)
        {
            const Individual &p1 = roulette_select_one(population, roulette_dist, rng);

            Individual &child = new_population[i];

            if (prob_dist(rng) < GeneticConfig::CROSSOVER_PROBABILITY)
            {
                roulette_dist_no_p1 = std::uniform_real_distribution<double>(0.0, fitness_sum - p1.fitness);
                const Individual &p2 = roulette_select_second_no_rep(population, p1, roulette_dist_no_p1, rng);

                crossover(p1.params, p2.params, child.params);
            }
            else
            {
                child = p1;
            }

            mutate(child.params, rng, ndist, prob_dist, esyn_mut_factor, use_ifast, use_islow);
        }

        std::swap(population_ptr, new_population_ptr);

        calc_fitnesses<Integrator, NeuronType, POP>(
            synapsis, model_neur, *population_ptr, scaled_result,
            living_const_sig_fitness_vals, search_phase, buffers,
            reset_state_neur, get_v_neur, ELITES, stabilization_points, avg_smooth_points,
            use_ifast, use_islow);
    }

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
