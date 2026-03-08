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
    static constexpr size_t POPULATION_SIZE = 30;
    static constexpr size_t NUM_GENERATIONS = 30;
    static constexpr size_t NUM_ELITES = 3;

    static constexpr double OBSERVATION_TIME_DIVISOR = 3.0;

    static constexpr double CROSSOVER_PROBABILITY = 0.9;
    static constexpr double MUTATION_PROBABILITY = 0.1;

    static constexpr double ETA = 0.2;

    static constexpr double ESYN_PHASE = -15.0;
    static constexpr double ESYN_ANTIPHASE = 15.0;

    static constexpr double SFAST_MIN = 1.25;
    static constexpr double SFAST_MAX = 4.0;

    static constexpr double SSLOW_MIN = 0.5;
    static constexpr double SSLOW_MAX = 1.0;

    static constexpr double K1_MIN = 0.0000001;
    static constexpr double K1_MAX = 1.0;

    static constexpr double K2_MIN = 0.0000001;
    static constexpr double K2_MAX = 1.0;

    static constexpr double GFAST_MIN = 0.001;
    static constexpr double GFAST_MAX = 3.0;

    static constexpr double GSLOW_MIN = 0.001;
    static constexpr double GSLOW_MAX = 5.0;

    static constexpr double VFAST_MIN = HindmarshRose::MIN;
    static constexpr double VFAST_MAX = HindmarshRose::MAX;

    static constexpr double VSLOW_PHASE_MIN = (HindmarshRose::MAX - HindmarshRose::MIN) * (3.0 / 4.0) + HindmarshRose::MIN;
    static constexpr double VSLOW_PHASE_MAX = HindmarshRose::MAX;

    static constexpr double VSLOW_ANTIPHASE_MIN = HindmarshRose::MIN;
    static constexpr double VSLOW_ANTIPHASE_MAX = (HindmarshRose::MAX - HindmarshRose::MIN) * (1.0 / 4.0) + HindmarshRose::MIN;

    static constexpr double SFAST_MUT_FACTOR = ETA * (SFAST_MAX - SFAST_MIN);
    static constexpr double VFAST_MUT_FACTOR = ETA * (VFAST_MAX - VFAST_MIN);
    static constexpr double VSLOW_PHASE_MUT_FACTOR = ETA * (VSLOW_PHASE_MAX - VSLOW_PHASE_MIN);
    static constexpr double VSLOW_ANTIPHASE_MUT_FACTOR = ETA * (VSLOW_ANTIPHASE_MAX - VSLOW_ANTIPHASE_MIN);
    static constexpr double K1_MUT_FACTOR = ETA * (K1_MAX - K1_MIN);
    static constexpr double K2_MUT_FACTOR = ETA * (K2_MAX - K2_MIN);
    static constexpr double SSLOW_MUT_FACTOR = ETA * (SSLOW_MAX - SSLOW_MIN);
    static constexpr double GFAST_MUT_FACTOR = ETA * (GFAST_MAX - GFAST_MIN);
    static constexpr double GSLOW_MUT_FACTOR = ETA * (GSLOW_MAX - GSLOW_MIN);
}

namespace GeneticConstants
{
    static const double LOG_K1_MIN = std::log(GeneticConfig::K1_MIN);
    static const double LOG_K1_MAX = std::log(GeneticConfig::K1_MAX);
    static const double LOG_K2_MIN = std::log(GeneticConfig::K2_MIN);
    static const double LOG_K2_MAX = std::log(GeneticConfig::K2_MAX);
    static const double LOG_GFAST_MIN = std::log(GeneticConfig::GFAST_MIN);
    static const double LOG_GFAST_MAX = std::log(GeneticConfig::GFAST_MAX);
    static const double LOG_GSLOW_MIN = std::log(GeneticConfig::GSLOW_MIN);
    static const double LOG_GSLOW_MAX = std::log(GeneticConfig::GSLOW_MAX);

    static const double LOG_K1_RANGE = LOG_K1_MAX - LOG_K1_MIN;
    static const double LOG_K2_RANGE = LOG_K2_MAX - LOG_K2_MIN;
    static const double LOG_GFAST_RANGE = LOG_GFAST_MAX - LOG_GFAST_MIN;
    static const double LOG_GSLOW_RANGE = LOG_GSLOW_MAX - LOG_GSLOW_MIN;
}

static inline double bounce_clamp(double value, double min, double max)
{
    while (value < min || value > max)
    {
        if (value < min)
            value = min + (min - value);
        if (value > max)
            value = max - (value - max);
    }
    return value;
}

static inline bool fitness_descending(const Individual &a, const Individual &b)
{
    return a.fitness > b.fitness;
}

template <size_t POP_SIZE>
static inline std::array<Individual, POP_SIZE> initialize_population(std::mt19937 &rng,
                                                                     bool use_ifast,
                                                                     bool use_islow,
                                                                     bool search_phase)
{
    std::uniform_real_distribution<double> dist_sfast;
    std::uniform_real_distribution<double> dist_vfast;
    std::uniform_real_distribution<double> dist_vslow;
    std::uniform_real_distribution<double> dist_k1;
    std::uniform_real_distribution<double> dist_k2;
    std::uniform_real_distribution<double> dist_sslow;
    std::uniform_real_distribution<double> dist_gfast;
    std::uniform_real_distribution<double> dist_gslow;

    if (use_ifast)
    {
        dist_sfast = std::uniform_real_distribution<double>(GeneticConfig::SFAST_MIN, GeneticConfig::SFAST_MAX);
        dist_vfast = std::uniform_real_distribution<double>(GeneticConfig::VFAST_MIN, GeneticConfig::VFAST_MAX);
        dist_gfast = std::uniform_real_distribution<double>(GeneticConstants::LOG_GFAST_MIN, GeneticConstants::LOG_GFAST_MAX);
    }

    if (use_islow)
    {
        dist_vslow = search_phase
                         ? std::uniform_real_distribution<double>(GeneticConfig::VSLOW_PHASE_MIN, GeneticConfig::VSLOW_PHASE_MAX)
                         : std::uniform_real_distribution<double>(GeneticConfig::VSLOW_ANTIPHASE_MIN, GeneticConfig::VSLOW_ANTIPHASE_MAX);
        dist_k1 = std::uniform_real_distribution<double>(GeneticConstants::LOG_K1_MIN, GeneticConstants::LOG_K1_MAX);
        dist_k2 = std::uniform_real_distribution<double>(GeneticConstants::LOG_K2_MIN, GeneticConstants::LOG_K2_MAX);
        dist_sslow = std::uniform_real_distribution<double>(GeneticConfig::SSLOW_MIN, GeneticConfig::SSLOW_MAX);
        dist_gslow = std::uniform_real_distribution<double>(GeneticConstants::LOG_GSLOW_MIN, GeneticConstants::LOG_GSLOW_MAX);
    }

    std::array<Individual, POP_SIZE> population;

    for (Individual &ind : population)
    {
        ChemicalSynapsisVariationParams &p = ind.params;
        if (use_ifast)
        {
            p.gfast = std::exp(dist_gfast(rng));
            p.sfast = dist_sfast(rng);
            p.Vfast = dist_vfast(rng);
        }

        if (use_islow)
        {
            p.gslow = std::exp(dist_gslow(rng));
            p.Vslow = dist_vslow(rng);
            p.k1 = std::exp(dist_k1(rng));
            p.k2 = std::exp(dist_k2(rng));
            p.sslow = dist_sslow(rng);
        }
    }

    return population;
}

static inline void crossover(const ChemicalSynapsisVariationParams &a, const ChemicalSynapsisVariationParams &b, ChemicalSynapsisVariationParams &result)
{
    result.gfast = (a.gfast + b.gfast) * 0.5;
    result.gslow = (a.gslow + b.gslow) * 0.5;
    result.sfast = (a.sfast + b.sfast) * 0.5;
    result.Vfast = (a.Vfast + b.Vfast) * 0.5;
    result.Vslow = (a.Vslow + b.Vslow) * 0.5;
    result.k1 = (a.k1 + b.k1) * 0.5;
    result.k2 = (a.k2 + b.k2) * 0.5;
    result.sslow = (a.sslow + b.sslow) * 0.5;
}

static inline void mutate(ChemicalSynapsisVariationParams &p, std::mt19937 &rng, std::normal_distribution<double> &ndist, std::uniform_real_distribution<double> &prob_dist, bool use_ifast, bool use_islow, bool search_phase)
{
    constexpr double MUTATION_PROBABILITY = GeneticConfig::MUTATION_PROBABILITY;
    constexpr double ETA = GeneticConfig::ETA;

    if (use_ifast)
    {
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &gfast = p.gfast;
            double logval = std::log(gfast);
            logval += ndist(rng) * (ETA * GeneticConstants::LOG_GFAST_RANGE);
            logval = bounce_clamp(logval, GeneticConstants::LOG_GFAST_MIN, GeneticConstants::LOG_GFAST_MAX);
            gfast = std::exp(logval);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &sfast = p.sfast;
            sfast += ndist(rng) * GeneticConfig::SFAST_MUT_FACTOR;
            sfast = bounce_clamp(sfast, GeneticConfig::SFAST_MIN, GeneticConfig::SFAST_MAX);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &Vfast = p.Vfast;
            Vfast += ndist(rng) * GeneticConfig::VFAST_MUT_FACTOR;
            Vfast = bounce_clamp(Vfast, GeneticConfig::VFAST_MIN, GeneticConfig::VFAST_MAX);
        }
    }

    if (use_islow)
    {
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &gslow = p.gslow;
            double logval = std::log(gslow);
            logval += ndist(rng) * (ETA * GeneticConstants::LOG_GSLOW_RANGE);
            logval = bounce_clamp(logval, GeneticConstants::LOG_GSLOW_MIN, GeneticConstants::LOG_GSLOW_MAX);
            gslow = std::exp(logval);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &Vslow = p.Vslow;
            double vslow_min, vslow_max, vslow_mut;
            if (search_phase)
            {
                vslow_min = GeneticConfig::VSLOW_PHASE_MIN;
                vslow_max = GeneticConfig::VSLOW_PHASE_MAX;
                vslow_mut = GeneticConfig::VSLOW_PHASE_MUT_FACTOR;
            }
            else
            {
                vslow_min = GeneticConfig::VSLOW_ANTIPHASE_MIN;
                vslow_max = GeneticConfig::VSLOW_ANTIPHASE_MAX;
                vslow_mut = GeneticConfig::VSLOW_ANTIPHASE_MUT_FACTOR;
            }
            Vslow += ndist(rng) * vslow_mut;
            Vslow = bounce_clamp(Vslow, vslow_min, vslow_max);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &k1 = p.k1;
            double logval = std::log(k1);
            logval += ndist(rng) * (ETA * GeneticConstants::LOG_K1_RANGE);
            logval = bounce_clamp(logval, GeneticConstants::LOG_K1_MIN, GeneticConstants::LOG_K1_MAX);
            k1 = std::exp(logval);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &k2 = p.k2;
            double logval = std::log(k2);
            logval += ndist(rng) * (ETA * GeneticConstants::LOG_K2_RANGE);
            logval = bounce_clamp(logval, GeneticConstants::LOG_K2_MIN, GeneticConstants::LOG_K2_MAX);
            k2 = std::exp(logval);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &sslow = p.sslow;
            sslow += ndist(rng) * GeneticConfig::SSLOW_MUT_FACTOR;
            sslow = bounce_clamp(sslow, GeneticConfig::SSLOW_MIN, GeneticConfig::SSLOW_MAX);
        }
    }
}

template <size_t N>
static inline const Individual &roulette_select_one(const std::array<Individual, N> &population,
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
static inline const Individual &roulette_select_second_no_rep(const std::array<Individual, N> &population,
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

    const bool use_ifast = (syn_component != SynComponent::ISLOW);
    const bool use_islow = (syn_component != SynComponent::IFAST);

    const double used_esyn = search_phase ? GeneticConfig::ESYN_PHASE : GeneticConfig::ESYN_ANTIPHASE;

    typename ChemicalSynapsisType::ConstructorArgs syn_args{};
    syn_args.params[ChemicalSynapsisType::Esyn] = used_esyn;
    if (use_islow && !use_ifast)
    {
        syn_args.params[ChemicalSynapsisType::gfast] = 0.0;
        syn_args.params[ChemicalSynapsisType::sfast] = 0.0;
        syn_args.params[ChemicalSynapsisType::Vfast] = 0.0;
    }
    else if (use_ifast && !use_islow)
    {
        syn_args.params[ChemicalSynapsisType::gslow] = 0.0;
        syn_args.params[ChemicalSynapsisType::Vslow] = 0.0;
        syn_args.params[ChemicalSynapsisType::k1] = 0.0;
        syn_args.params[ChemicalSynapsisType::k2] = 0.0;
        syn_args.params[ChemicalSynapsisType::sslow] = 0.0;
    }

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

    const double pts_burst_real = scaled_result.pts_burst_real;
    const size_t avg_smooth_points = std::max(FitnessConstants::MIN_AVG_SMOOTH_POINTS, static_cast<size_t>(pts_burst_real / FitnessConfig::AVG_SMOOTH_POINTS_BURST_DIVISOR));

    const size_t stabilization_points = std::max(static_cast<size_t>(stabilization_time / csv_step), avg_smooth_points);
    const size_t use_vpre_sig_size = scaled_result.sig.size() - stabilization_points;

    const size_t smoothing_size = avg_smooth_points + use_vpre_sig_size;

    SigBuffers buffers;
    buffers.vpost_sig.resize(smoothing_size);
    if (use_ifast && use_islow)
        buffers.i_sig.resize(use_vpre_sig_size);
    if (use_ifast)
        buffers.ifast_sig.resize(use_vpre_sig_size);
    if (use_islow)
        buffers.islow_sig.resize(use_vpre_sig_size);
    buffers.smoothed_vpost_sig.resize(use_vpre_sig_size);

    const ConstantSigFitnessVals const_vpre_sig_fitness_vals = calc_const_sig_fitness_vals(
        scaled_result.sig.slice(stabilization_points - avg_smooth_points, smoothing_size),
        model_min, model_max, search_phase,
        avg_smooth_points, pts_burst_real, use_ifast, use_islow);

    constexpr size_t POP = GeneticConfig::POPULATION_SIZE;
    constexpr size_t ELITES = GeneticConfig::NUM_ELITES;

    std::random_device rd;
    std::mt19937 rng(rd());

    std::array<Individual, POP> buffer1 = initialize_population<POP>(rng, use_ifast, use_islow, search_phase);
    std::array<Individual, POP> buffer2;

    std::array<Individual, POP> *population_ptr = &buffer1;
    std::array<Individual, POP> *new_population_ptr = &buffer2;

    calc_fitnesses<Integrator, NeuronType, POP>(
        synapsis, model_neur, *population_ptr, scaled_result,
        const_vpre_sig_fitness_vals, search_phase, buffers,
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

        std::copy_n(population.begin(), ELITES, new_population.begin());

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

            mutate(child.params, rng, ndist, prob_dist, use_ifast, use_islow, search_phase);
        }

        std::swap(population_ptr, new_population_ptr);

        calc_fitnesses<Integrator, NeuronType, POP>(
            synapsis, model_neur, *population_ptr, scaled_result,
            const_vpre_sig_fitness_vals, search_phase, buffers,
            reset_state_neur, get_v_neur, ELITES, stabilization_points, avg_smooth_points,
            use_ifast, use_islow);
    }

    std::array<Individual, POP> &population = *population_ptr;
    const Individual &best = *std::min_element(population.begin(), population.end(), fitness_descending);

    if (verbose)
    {
        std::cout << std::endl;
        const ChemicalSynapsisVariationParams &params = best.params;
        std::cout << "Fitness: " << best.fitness << std::endl;
        std::cout << "gfast_values = [" << params.gfast << "]" << std::endl;
        std::cout << "gslow_values = [" << params.gslow << "]" << std::endl;
        std::cout << "Esyn_values = [" << used_esyn << "]" << std::endl;
        std::cout << "Vfast_values = [" << params.Vfast << "]" << std::endl;
        std::cout << "Vslow_values = [" << params.Vslow << "]" << std::endl;
        std::cout << "sfast_values = [" << params.sfast << "]" << std::endl;
        std::cout << "sslow_values = [" << params.sslow << "]" << std::endl;
        std::cout << "k1_values = [" << params.k1 << "]" << std::endl;
        std::cout << "k2_values = [" << params.k2 << "]" << std::endl;
    }

    return best;
}
