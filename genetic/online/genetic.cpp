#include "bidirectional_chemical_synapse_genetic.h"
#include <algorithm>
#include <cmath>
#include "utils.hpp"

namespace GeneticPrivateConfig
{
    static constexpr size_t NUM_ELITES = 1;

    static constexpr double CROSSOVER_PROBABILITY = 0.9;
    static constexpr double MUTATION_PROBABILITY = 0.1;
}

static double bounce_clamp(double value, double min, double max)
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

static inline const Individual &roulette_select_one(const std::vector<Individual> &population,
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

static inline const Individual &roulette_select_second_no_rep(const std::vector<Individual> &population,
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

static void crossover_one_direction(const ChemicalSynapseVariationParams &a,
                                    const ChemicalSynapseVariationParams &b,
                                    ChemicalSynapseVariationParams &result,
                                    unsigned int use_i_fast,
                                    unsigned int use_i_slow)
{
    if (use_i_fast)
    {
        result.g_fast = (a.g_fast + b.g_fast) * 0.5;
        result.s_fast = (a.s_fast + b.s_fast) * 0.5;
        result.v_fast = (a.v_fast + b.v_fast) * 0.5;
    }

    if (use_i_slow)
    {
        result.g_slow = (a.g_slow + b.g_slow) * 0.5;
        result.v_slow = (a.v_slow + b.v_slow) * 0.5;
        result.k1 = (a.k1 + b.k1) * 0.5;
        result.k2 = (a.k2 + b.k2) * 0.5;
        result.s_slow = (a.s_slow + b.s_slow) * 0.5;
    }
}

static void mutate_one_direction(ChemicalSynapseVariationParams &p,
                                 std::mt19937 &rng,
                                 std::normal_distribution<double> &ndist,
                                 std::uniform_real_distribution<double> &prob_dist,
                                 unsigned int use_i_fast,
                                 unsigned int use_i_slow,
                                 const GeneticRanges &ranges)
{
    constexpr double MUTATION_PROBABILITY = GeneticPrivateConfig::MUTATION_PROBABILITY;

    if (use_i_fast)
    {
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &g_fast = p.g_fast;
            double logval = std::log(g_fast);
            logval += ndist(rng) * ranges.log_g_fast.mut_factor;
            logval = bounce_clamp(logval, ranges.log_g_fast.min, ranges.log_g_fast.max);
            g_fast = std::exp(logval);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &s_fast = p.s_fast;
            s_fast += ndist(rng) * ranges.s_fast.mut_factor;
            s_fast = bounce_clamp(s_fast, ranges.s_fast.min, ranges.s_fast.max);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &v_fast = p.v_fast;
            v_fast += ndist(rng) * ranges.v_fast.mut_factor;
            v_fast = bounce_clamp(v_fast, ranges.v_fast.min, ranges.v_fast.max);
        }
    }

    if (use_i_slow)
    {
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &g_slow = p.g_slow;
            double logval = std::log(g_slow);
            logval += ndist(rng) * ranges.log_g_slow.mut_factor;
            logval = bounce_clamp(logval, ranges.log_g_slow.min, ranges.log_g_slow.max);
            g_slow = std::exp(logval);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &v_slow = p.v_slow;
            v_slow += ndist(rng) * ranges.v_slow.mut_factor;
            v_slow = bounce_clamp(v_slow, ranges.v_slow.min, ranges.v_slow.max);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &k1 = p.k1;
            double logval = std::log(k1);
            logval += ndist(rng) * ranges.log_k1.mut_factor;
            logval = bounce_clamp(logval, ranges.log_k1.min, ranges.log_k1.max);
            k1 = std::exp(logval);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &k2 = p.k2;
            double logval = std::log(k2);
            logval += ndist(rng) * ranges.log_k2.mut_factor;
            logval = bounce_clamp(logval, ranges.log_k2.min, ranges.log_k2.max);
            k2 = std::exp(logval);
        }
        if (prob_dist(rng) < MUTATION_PROBABILITY)
        {
            double &s_slow = p.s_slow;
            s_slow += ndist(rng) * ranges.s_slow.mut_factor;
            s_slow = bounce_clamp(s_slow, ranges.s_slow.min, ranges.s_slow.max);
        }
    }
}

struct OneDirectionInitializeDistributions
{
    std::uniform_real_distribution<double> dist_s_fast;
    std::uniform_real_distribution<double> dist_v_fast;
    std::uniform_real_distribution<double> dist_log_g_fast;
    std::uniform_real_distribution<double> dist_s_slow;
    std::uniform_real_distribution<double> dist_v_slow;
    std::uniform_real_distribution<double> dist_log_k1;
    std::uniform_real_distribution<double> dist_log_k2;
    std::uniform_real_distribution<double> dist_log_g_slow;

    OneDirectionInitializeDistributions(unsigned int use_i_fast,
                                        unsigned int use_i_slow,
                                        const GeneticRanges &ranges)
    {
        if (use_i_fast)
        {
            dist_s_fast = std::uniform_real_distribution<double>(ranges.s_fast.min, ranges.s_fast.max);
            dist_v_fast = std::uniform_real_distribution<double>(ranges.v_fast.min, ranges.v_fast.max);
            dist_log_g_fast = std::uniform_real_distribution<double>(ranges.log_g_fast.min, ranges.log_g_fast.max);
        }

        if (use_i_slow)
        {
            dist_s_slow = std::uniform_real_distribution<double>(ranges.s_slow.min, ranges.s_slow.max);
            dist_v_slow = std::uniform_real_distribution<double>(ranges.v_slow.min, ranges.v_slow.max);
            dist_log_k1 = std::uniform_real_distribution<double>(ranges.log_k1.min, ranges.log_k1.max);
            dist_log_k2 = std::uniform_real_distribution<double>(ranges.log_k2.min, ranges.log_k2.max);
            dist_log_g_slow = std::uniform_real_distribution<double>(ranges.log_g_slow.min, ranges.log_g_slow.max);
        }
    }
};

static void initialize_individual_one_direction(ChemicalSynapseVariationParams &p,
                                                std::mt19937 &rng,
                                                unsigned int use_i_fast,
                                                unsigned int use_i_slow,
                                                OneDirectionInitializeDistributions &dists)
{
    if (use_i_fast)
    {
        p.g_fast = std::exp(dists.dist_log_g_fast(rng));
        p.s_fast = dists.dist_s_fast(rng);
        p.v_fast = dists.dist_v_fast(rng);
    }

    if (use_i_slow)
    {
        p.g_slow = std::exp(dists.dist_log_g_slow(rng));
        p.v_slow = dists.dist_v_slow(rng);
        p.k1 = std::exp(dists.dist_log_k1(rng));
        p.k2 = std::exp(dists.dist_log_k2(rng));
        p.s_slow = dists.dist_s_slow(rng);
    }
}

inline std::vector<Individual> BidirectionalChemicalSynapseGenetic::initialize_population(std::mt19937 &rng,
                                                                                          const GeneticRanges &ranges_12,
                                                                                          const GeneticRanges &ranges_21)
{
    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

    OneDirectionInitializeDistributions dists_12(use_i_fast_12, use_i_slow_12, ranges_12);
    OneDirectionInitializeDistributions dists_21(use_i_fast_21, use_i_slow_21, ranges_21);

    std::vector<Individual> population(population_size);

    for (Individual &ind : population)
    {
        if (use_syn_12)
            initialize_individual_one_direction(ind.variation_params_12, rng, use_i_fast_12, use_i_slow_12, dists_12);
        if (use_syn_21)
            initialize_individual_one_direction(ind.variation_params_21, rng, use_i_fast_21, use_i_slow_21, dists_21);
    }

    return population;
}

inline void BidirectionalChemicalSynapseGenetic::crossover_individual(const Individual &a, const Individual &b, Individual &result)
{
    if (use_i_fast_12 || use_i_slow_12)
        crossover_one_direction(a.variation_params_12, b.variation_params_12, result.variation_params_12,
                                use_i_fast_12, use_i_slow_12);

    if (use_i_fast_21 || use_i_slow_21)
        crossover_one_direction(a.variation_params_21, b.variation_params_21, result.variation_params_21,
                                use_i_fast_21, use_i_slow_21);
}

inline void BidirectionalChemicalSynapseGenetic::mutate_individual(Individual &ind,
                                                                   std::mt19937 &rng,
                                                                   std::normal_distribution<double> &ndist,
                                                                   std::uniform_real_distribution<double> &prob_dist,
                                                                   const GeneticRanges &ranges_12, const GeneticRanges &ranges_21)
{
    if (use_i_fast_12 || use_i_slow_12)
        mutate_one_direction(ind.variation_params_12, rng, ndist, prob_dist, use_i_fast_12, use_i_slow_12, ranges_12);
    if (use_i_fast_21 || use_i_slow_21)
        mutate_one_direction(ind.variation_params_21, rng, ndist, prob_dist, use_i_fast_21, use_i_slow_21, ranges_21);
}

void BidirectionalChemicalSynapseGenetic::NRT_genetic(double period_t)
{
    const double fs = 1.0 / period_t;                         // period está en ms, fs en KHz
    num_elements = static_cast<size_t>(evaluation_time * fs); // evaluation_time en ms y fs en KHz

    const size_t effective_pad = std::min(num_elements - 1,
                                          static_cast<size_t>(FitnessPublicConfig::PAD_LEN_FACTOR * fs / FitnessPublicConfig::FILTER_FC));
    univector<double> padded_buff(num_elements + (2 * effective_pad));

    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

    size_t curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);
    if (!wait_until_RT_read_idx_or_stop(curr_synapse_idx))
    {
        QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
        return;
    }
    size_t new_synapse_idx = 1 - curr_synapse_idx;

    GeneticRanges ranges_12;
    if (use_syn_12)
    {
        ranges_12.init(v_min_1,
                       v_max_1,
                       v_min_2,
                       v_max_2,
                       expected_i_min_12,
                       expected_i_max_12,
                       use_i_fast_12,
                       use_i_slow_12,
                       search_phase);
        params_12[new_synapse_idx].e_syn = ranges_12.e_syn;
        v_sig_1.resize(num_elements);
        if (use_i_fast_12)
            i_fast_sig_12.resize(num_elements);
        if (use_i_slow_12)
            i_slow_sig_12.resize(num_elements);
    }

    GeneticRanges ranges_21;
    if (use_syn_21)
    {
        ranges_21.init(v_min_2,
                       v_max_2,
                       v_min_1,
                       v_max_1,
                       expected_i_min_21,
                       expected_i_max_21,
                       use_i_fast_21,
                       use_i_slow_21,
                       search_phase);
        params_21[new_synapse_idx].e_syn = ranges_21.e_syn;
        v_sig_2.resize(num_elements);
        if (use_i_fast_21)
            i_fast_sig_21.resize(num_elements);
        if (use_i_slow_21)
            i_slow_sig_21.resize(num_elements);
    }

    synapse_idx.store(new_synapse_idx, std::memory_order_release);
    if (!wait_until_RT_read_idx_or_stop(new_synapse_idx))
    {
        QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
        return;
    }

    if (use_syn_12)
        params_12[curr_synapse_idx].e_syn = ranges_12.e_syn;
    if (use_syn_21)
        params_21[curr_synapse_idx].e_syn = ranges_21.e_syn;

    synapse_idx.store(curr_synapse_idx, std::memory_order_release);

    constexpr size_t ELITES = GeneticPrivateConfig::NUM_ELITES;

    std::random_device rd;
    std::mt19937 rng(rd());

    std::vector<Individual> population = initialize_population(rng, ranges_12, ranges_21);
    std::vector<Individual> new_population(population_size);

    QMetaObject::invokeMethod(this, "set_generations_completed", Qt::QueuedConnection,
                              Q_ARG(double, 0));
    QMetaObject::invokeMethod(this, "set_individuals_completed", Qt::QueuedConnection,
                              Q_ARG(double, 0));

    if (!calc_fitnesses(std::span<Individual>(population), fs, effective_pad, padded_buff))
    {
        QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
        return;
    }

    QMetaObject::invokeMethod(this, "set_generations_completed", Qt::QueuedConnection,
                              Q_ARG(double, 1));

    std::normal_distribution<double> ndist(0.0, 1.0);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
    std::uniform_real_distribution<double> roulette_dist;
    std::uniform_real_distribution<double> roulette_dist_no_p1;

    typename std::vector<Individual>::iterator pop_begin;

    for (size_t gen = 0; gen < num_generations; gen++)
    {
        QMetaObject::invokeMethod(this, "set_individuals_completed", Qt::QueuedConnection,
                                  Q_ARG(double, 0));

        if (stop_genetic.load(std::memory_order_relaxed))
        {
            QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
            return;
        }

        pop_begin = population.begin();
        std::nth_element(pop_begin, pop_begin + ELITES, population.end(), fitness_descending);

        double fitness_sum = 0.0;
        for (const Individual &ind : population)
            fitness_sum += ind.fitness;

        roulette_dist = std::uniform_real_distribution<double>(0.0, fitness_sum);

        std::copy_n(pop_begin, ELITES, new_population.begin());

        QMetaObject::invokeMethod(this, "set_individuals_completed", Qt::QueuedConnection,
                                  Q_ARG(double, ELITES));

        for (size_t i = ELITES; i < population_size; i++)
        {
            const Individual &p1 = roulette_select_one(population, roulette_dist, rng);

            Individual &child = new_population[i];

            if (prob_dist(rng) < GeneticPrivateConfig::CROSSOVER_PROBABILITY)
            {
                roulette_dist_no_p1 = std::uniform_real_distribution<double>(0.0, fitness_sum - p1.fitness);
                const Individual &p2 = roulette_select_second_no_rep(population, p1, roulette_dist_no_p1, rng);

                crossover_individual(p1, p2, child);
            }
            else
            {
                child = p1;
            }

            mutate_individual(child, rng, ndist, prob_dist, ranges_12, ranges_21);
        }

        std::swap(population, new_population);

        if (!calc_fitnesses(std::span<Individual>(population).subspan(ELITES), fs, effective_pad, padded_buff))
        {
            QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
            return;
        }

        QMetaObject::invokeMethod(this, "set_generations_completed", Qt::QueuedConnection,
                                  Q_ARG(double, gen + 1));
    }

    const Individual &best = *std::min_element(population.begin(), population.end(), fitness_descending);

    curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);
    if (!wait_until_RT_read_idx_or_stop(curr_synapse_idx))
    {
        QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
        return;
    }

    new_synapse_idx = 1 - curr_synapse_idx;
    if (use_syn_12)
    {
        apply_variation_params(params_12[new_synapse_idx],
                               best.variation_params_12,
                               use_i_fast_12,
                               use_i_slow_12);
    }
    if (use_syn_21)
    {
        apply_variation_params(params_21[new_synapse_idx],
                               best.variation_params_21,
                               use_i_fast_21,
                               use_i_slow_21);
    }
    synapse_idx.store(new_synapse_idx, std::memory_order_release);

    QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
}
