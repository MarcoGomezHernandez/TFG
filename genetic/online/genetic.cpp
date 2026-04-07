#include "bidirectional_chemical_synapse_genetic.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include "utils.hpp"

namespace GeneticPrivateConfig
{
    static constexpr size_t NUM_ELITES = 1;
    static constexpr unsigned int TOUR_K = 3;

    static constexpr double CROSSOVER_PROBABILITY = 0.9;
    static constexpr double INDIVIDUAL_MUTATION_PROBABILITY = 0.2;
    static constexpr double BLX_ALPHA = 0.1;
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

static void crossover_blx_alpha_param(double a,
                                      double b,
                                      const GeneticRanges::ParamRange &range,
                                      std::mt19937 &rng,
                                      double &result_1,
                                      double &result_2,
                                      bool has_second_individual)
{
    const double min = std::min(a, b);
    const double max = std::max(a, b);
    const double interval = max - min;

    constexpr double BLX_ALPHA = GeneticPrivateConfig::BLX_ALPHA;

    const double lower = min - (BLX_ALPHA * interval);
    const double upper = max + (BLX_ALPHA * interval);

    std::uniform_real_distribution<double> dist(lower, upper);
    result_1 = bounce_clamp(dist(rng), range.min, range.max);
    if (has_second_individual)
        result_2 = bounce_clamp(dist(rng), range.min, range.max);
}

static size_t num_mutation_params_one_direction(unsigned int use_i_fast,
                                                unsigned int use_i_slow)
{
    size_t num_params = 0;
    if (use_i_fast || use_i_slow)
    {
        num_params += 1;
        if (use_i_fast)
            num_params += 3;
        if (use_i_slow)
            num_params += 5;
    }

    return num_params;
}

static void crossover_one_direction(const IndividualChemicalSynapseParams &a,
                                    const IndividualChemicalSynapseParams &b,
                                    IndividualChemicalSynapseParams &result_1,
                                    IndividualChemicalSynapseParams &result_2,
                                    bool has_second_individual,
                                    unsigned int use_i_fast,
                                    unsigned int use_i_slow,
                                    const GeneticRanges &ranges,
                                    std::mt19937 &rng)
{
    if (use_i_fast || use_i_slow)
    {
        crossover_blx_alpha_param(a.e_syn,
                                  b.e_syn,
                                  ranges.e_syn,
                                  rng,
                                  result_1.e_syn,
                                  result_2.e_syn,
                                  has_second_individual);

        if (use_i_fast)
        {
            crossover_blx_alpha_param(a.log_g_fast,
                                      b.log_g_fast,
                                      ranges.log_g_fast,
                                      rng,
                                      result_1.log_g_fast,
                                      result_2.log_g_fast,
                                      has_second_individual);
            crossover_blx_alpha_param(a.s_fast,
                                      b.s_fast,
                                      ranges.s_fast,
                                      rng,
                                      result_1.s_fast,
                                      result_2.s_fast,
                                      has_second_individual);
            crossover_blx_alpha_param(a.v_fast,
                                      b.v_fast,
                                      ranges.v_fast,
                                      rng,
                                      result_1.v_fast,
                                      result_2.v_fast,
                                      has_second_individual);
        }

        if (use_i_slow)
        {
            crossover_blx_alpha_param(a.log_g_slow,
                                      b.log_g_slow,
                                      ranges.log_g_slow,
                                      rng,
                                      result_1.log_g_slow,
                                      result_2.log_g_slow,
                                      has_second_individual);
            crossover_blx_alpha_param(a.v_slow,
                                      b.v_slow,
                                      ranges.v_slow,
                                      rng,
                                      result_1.v_slow,
                                      result_2.v_slow,
                                      has_second_individual);
            crossover_blx_alpha_param(a.log_k1,
                                      b.log_k1,
                                      ranges.log_k1,
                                      rng,
                                      result_1.log_k1,
                                      result_2.log_k1,
                                      has_second_individual);
            crossover_blx_alpha_param(a.log_k2,
                                      b.log_k2,
                                      ranges.log_k2,
                                      rng,
                                      result_1.log_k2,
                                      result_2.log_k2,
                                      has_second_individual);
            crossover_blx_alpha_param(a.s_slow,
                                      b.s_slow,
                                      ranges.s_slow,
                                      rng,
                                      result_1.s_slow,
                                      result_2.s_slow,
                                      has_second_individual);
        }
    }
}

static void mutate_one_direction(IndividualChemicalSynapseParams &p,
                                 std::mt19937 &rng,
                                 std::normal_distribution<double> &ndist,
                                 std::uniform_real_distribution<double> &prob_dist,
                                 unsigned int use_i_fast,
                                 unsigned int use_i_slow,
                                 double mutation_probability_per_gene,
                                 const GeneticRanges &ranges)
{
    if (use_i_fast || use_i_slow)
    {
        if (prob_dist(rng) < mutation_probability_per_gene)
        {
            const GeneticRanges::ParamRange &e_syn_range = ranges.e_syn;
            double &e_syn = p.e_syn;
            e_syn += ndist(rng) * e_syn_range.mut_factor;
            e_syn = bounce_clamp(e_syn, e_syn_range.min, e_syn_range.max);
        }

        if (use_i_fast)
        {
            if (prob_dist(rng) < mutation_probability_per_gene)
            {
                const GeneticRanges::ParamRange &log_g_fast_range = ranges.log_g_fast;
                double &log_g_fast = p.log_g_fast;
                log_g_fast += ndist(rng) * log_g_fast_range.mut_factor;
                log_g_fast = bounce_clamp(log_g_fast, log_g_fast_range.min, log_g_fast_range.max);
            }
            if (prob_dist(rng) < mutation_probability_per_gene)
            {
                const GeneticRanges::ParamRange &s_fast_range = ranges.s_fast;
                double &s_fast = p.s_fast;
                s_fast += ndist(rng) * s_fast_range.mut_factor;
                s_fast = bounce_clamp(s_fast, s_fast_range.min, s_fast_range.max);
            }
            if (prob_dist(rng) < mutation_probability_per_gene)
            {
                const GeneticRanges::ParamRange &v_fast_range = ranges.v_fast;
                double &v_fast = p.v_fast;
                v_fast += ndist(rng) * v_fast_range.mut_factor;
                v_fast = bounce_clamp(v_fast, v_fast_range.min, v_fast_range.max);
            }
        }

        if (use_i_slow)
        {
            if (prob_dist(rng) < mutation_probability_per_gene)
            {
                const GeneticRanges::ParamRange &log_g_slow_range = ranges.log_g_slow;
                double &log_g_slow = p.log_g_slow;
                log_g_slow += ndist(rng) * log_g_slow_range.mut_factor;
                log_g_slow = bounce_clamp(log_g_slow, log_g_slow_range.min, log_g_slow_range.max);
            }
            if (prob_dist(rng) < mutation_probability_per_gene)
            {
                const GeneticRanges::ParamRange &v_slow_range = ranges.v_slow;
                double &v_slow = p.v_slow;
                v_slow += ndist(rng) * v_slow_range.mut_factor;
                v_slow = bounce_clamp(v_slow, v_slow_range.min, v_slow_range.max);
            }
            if (prob_dist(rng) < mutation_probability_per_gene)
            {
                const GeneticRanges::ParamRange &log_k1_range = ranges.log_k1;
                double &log_k1 = p.log_k1;
                log_k1 += ndist(rng) * log_k1_range.mut_factor;
                log_k1 = bounce_clamp(log_k1, log_k1_range.min, log_k1_range.max);
            }
            if (prob_dist(rng) < mutation_probability_per_gene)
            {
                const GeneticRanges::ParamRange &log_k2_range = ranges.log_k2;
                double &log_k2 = p.log_k2;
                log_k2 += ndist(rng) * log_k2_range.mut_factor;
                log_k2 = bounce_clamp(log_k2, log_k2_range.min, log_k2_range.max);
            }
            if (prob_dist(rng) < mutation_probability_per_gene)
            {
                const GeneticRanges::ParamRange &s_slow_range = ranges.s_slow;
                double &s_slow = p.s_slow;
                s_slow += ndist(rng) * s_slow_range.mut_factor;
                s_slow = bounce_clamp(s_slow, s_slow_range.min, s_slow_range.max);
            }
        }
    }
}

struct OneDirectionInitializeDistributions
{
    std::uniform_real_distribution<double> dist_e_syn;
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
        if (use_i_fast || use_i_slow)
        {
            const GeneticRanges::ParamRange &e_syn_range = ranges.e_syn;
            dist_e_syn = std::uniform_real_distribution<double>(e_syn_range.min, e_syn_range.max);

            if (use_i_fast)
            {
                const GeneticRanges::ParamRange &s_fast_range = ranges.s_fast;
                const GeneticRanges::ParamRange &v_fast_range = ranges.v_fast;
                const GeneticRanges::ParamRange &log_g_fast_range = ranges.log_g_fast;

                dist_s_fast = std::uniform_real_distribution<double>(s_fast_range.min, s_fast_range.max);
                dist_v_fast = std::uniform_real_distribution<double>(v_fast_range.min, v_fast_range.max);
                dist_log_g_fast = std::uniform_real_distribution<double>(log_g_fast_range.min, log_g_fast_range.max);
            }

            if (use_i_slow)
            {
                const GeneticRanges::ParamRange &s_slow_range = ranges.s_slow;
                const GeneticRanges::ParamRange &v_slow_range = ranges.v_slow;
                const GeneticRanges::ParamRange &log_k1_range = ranges.log_k1;
                const GeneticRanges::ParamRange &log_k2_range = ranges.log_k2;
                const GeneticRanges::ParamRange &log_g_slow_range = ranges.log_g_slow;

                dist_s_slow = std::uniform_real_distribution<double>(s_slow_range.min, s_slow_range.max);
                dist_v_slow = std::uniform_real_distribution<double>(v_slow_range.min, v_slow_range.max);
                dist_log_k1 = std::uniform_real_distribution<double>(log_k1_range.min, log_k1_range.max);
                dist_log_k2 = std::uniform_real_distribution<double>(log_k2_range.min, log_k2_range.max);
                dist_log_g_slow = std::uniform_real_distribution<double>(log_g_slow_range.min, log_g_slow_range.max);
            }
        }
    }
};

static void initialize_individual_one_direction(IndividualChemicalSynapseParams &p,
                                                std::mt19937 &rng,
                                                unsigned int use_i_fast,
                                                unsigned int use_i_slow,
                                                OneDirectionInitializeDistributions &dists)
{
    if (use_i_fast || use_i_slow)
    {
        p.e_syn = dists.dist_e_syn(rng);

        if (use_i_fast)
        {
            p.log_g_fast = dists.dist_log_g_fast(rng);
            p.s_fast = dists.dist_s_fast(rng);
            p.v_fast = dists.dist_v_fast(rng);
        }

        if (use_i_slow)
        {
            p.log_g_slow = dists.dist_log_g_slow(rng);
            p.v_slow = dists.dist_v_slow(rng);
            p.log_k1 = dists.dist_log_k1(rng);
            p.log_k2 = dists.dist_log_k2(rng);
            p.s_slow = dists.dist_s_slow(rng);
        }
    }
}

inline std::vector<Individual> BidirectionalChemicalSynapseGenetic::initialize_population(std::mt19937 &rng,
                                                                                          const GeneticRanges &ranges_12,
                                                                                          const GeneticRanges &ranges_21)
{
    OneDirectionInitializeDistributions dists_12(use_i_fast_12, use_i_slow_12, ranges_12);
    OneDirectionInitializeDistributions dists_21(use_i_fast_21, use_i_slow_21, ranges_21);

    std::vector<Individual> population(population_size);

    for (Individual &ind : population)
    {
        initialize_individual_one_direction(ind.params_12, rng, use_i_fast_12, use_i_slow_12, dists_12);
        initialize_individual_one_direction(ind.params_21, rng, use_i_fast_21, use_i_slow_21, dists_21);
    }

    return population;
}

inline void BidirectionalChemicalSynapseGenetic::crossover_individual(const Individual &a,
                                                                      const Individual &b,
                                                                      Individual &result_1,
                                                                      Individual &result_2,
                                                                      bool has_second_individual,
                                                                      const GeneticRanges &ranges_12,
                                                                      const GeneticRanges &ranges_21,
                                                                      std::mt19937 &rng)
{
    crossover_one_direction(a.params_12,
                            b.params_12,
                            result_1.params_12,
                            result_2.params_12,
                            has_second_individual,
                            use_i_fast_12,
                            use_i_slow_12,
                            ranges_12,
                            rng);
    crossover_one_direction(a.params_21,
                            b.params_21,
                            result_1.params_21,
                            result_2.params_21,
                            has_second_individual,
                            use_i_fast_21,
                            use_i_slow_21,
                            ranges_21,
                            rng);
}

inline void BidirectionalChemicalSynapseGenetic::mutate_individual(Individual &ind,
                                                                   std::mt19937 &rng,
                                                                   std::normal_distribution<double> &ndist,
                                                                   std::uniform_real_distribution<double> &prob_dist,
                                                                   double mutation_probability_per_gene,
                                                                   const GeneticRanges &ranges_12, const GeneticRanges &ranges_21)
{
    mutate_one_direction(ind.params_12, rng, ndist, prob_dist, use_i_fast_12, use_i_slow_12,
                         mutation_probability_per_gene, ranges_12);
    mutate_one_direction(ind.params_21, rng, ndist, prob_dist, use_i_fast_21, use_i_slow_21,
                         mutation_probability_per_gene, ranges_21);
}

void BidirectionalChemicalSynapseGenetic::NRT_genetic(double period_t)
{
    const double fs = 1.0 / period_t;                         // period está en ms, fs en KHz
    num_elements = static_cast<size_t>(evaluation_time * fs); // evaluation_time en ms y fs en KHz
    const size_t pad_len_factor_fs = static_cast<size_t>(FitnessPublicConfig::PAD_LEN_FACTOR * fs);

    constexpr unsigned int TOUR_K = GeneticPrivateConfig::TOUR_K;

    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

    size_t effective_pad_12 = 0;
    size_t padded_buff_size_12 = 0;
    if (use_syn_12)
    {
        effective_pad_12 = std::min(num_elements - 1,
                                    static_cast<size_t>(pad_len_factor_fs / fc_1));
        padded_buff_size_12 = num_elements + (2 * effective_pad_12);
    }

    size_t effective_pad_21 = 0;
    size_t padded_buff_size_21 = 0;
    if (use_syn_21)
    {
        effective_pad_21 = std::min(num_elements - 1,
                                    static_cast<size_t>(pad_len_factor_fs / fc_2));
        padded_buff_size_21 = num_elements + (2 * effective_pad_21);
    }

    FitnessPadBuffers pad_buffers(padded_buff_size_12,
                                  padded_buff_size_21);

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
                       search_phase,
                       fc_1);
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
                       search_phase,
                       fc_2);
        v_sig_2.resize(num_elements);
        if (use_i_fast_21)
            i_fast_sig_21.resize(num_elements);
        if (use_i_slow_21)
            i_slow_sig_21.resize(num_elements);
    }

    constexpr size_t ELITES = GeneticPrivateConfig::NUM_ELITES;

    std::random_device rd;
    std::mt19937 rng(rd());

    std::vector<Individual> population = initialize_population(rng, ranges_12, ranges_21);
    std::vector<Individual> new_population(population_size);

    const unsigned int tour_k_p1 = std::min(TOUR_K, population_size);
    const unsigned int tour_k_p2 = population_size > 1 ? std::min(TOUR_K, population_size - 1) : 1u;

    std::vector<size_t> tour_idxs(population_size);
    std::vector<size_t>::iterator tour_idxs_begin = tour_idxs.begin();
    std::vector<size_t>::iterator tour_idxs_end = tour_idxs.end();
    std::iota(tour_idxs_begin, tour_idxs_end, static_cast<size_t>(0));

    std::vector<size_t> tour_sampled_idxs_p1(tour_k_p1);
    std::vector<size_t> tour_sampled_idxs_p2(tour_k_p2);

    std::vector<size_t>::iterator tour_sampled_idxs_p1_begin = tour_sampled_idxs_p1.begin();
    std::vector<size_t>::iterator tour_sampled_idxs_p2_begin = tour_sampled_idxs_p2.begin();

    std::function<bool(size_t, size_t)> fitness_less_by_population_index = [&population](size_t a, size_t b)
    {
        return population[a].fitness < population[b].fitness;
    };

    QMetaObject::invokeMethod(this, "set_generations_completed", Qt::QueuedConnection,
                              Q_ARG(double, 0));
    QMetaObject::invokeMethod(this, "set_individuals_completed", Qt::QueuedConnection,
                              Q_ARG(double, 0));

    if (!calc_fitnesses(std::span<Individual>(population),
                        fs,
                        effective_pad_12,
                        effective_pad_21,
                        pad_buffers))
    {
        QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
        return;
    }

    QMetaObject::invokeMethod(this, "set_generations_completed", Qt::QueuedConnection,
                              Q_ARG(double, 1));

    std::normal_distribution<double> ndist(0.0, 1.0);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

    const size_t num_params_12 = num_mutation_params_one_direction(use_i_fast_12, use_i_slow_12);
    const size_t num_params_21 = num_mutation_params_one_direction(use_i_fast_21, use_i_slow_21);
    const size_t num_params_total = num_params_12 + num_params_21;
    const double mutation_probability_per_gene = num_params_total > 0 ? (GeneticPrivateConfig::INDIVIDUAL_MUTATION_PROBABILITY / static_cast<double>(num_params_total)) : 0.0;

    Individual odd_aux_individual{};

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

        std::copy_n(pop_begin, ELITES, new_population.begin());

        QMetaObject::invokeMethod(this, "set_individuals_completed", Qt::QueuedConnection,
                                  Q_ARG(double, ELITES));

        for (size_t i = ELITES; i < population_size; i += 2)
        {
            std::sample(tour_idxs_begin,
                        tour_idxs_end,
                        tour_sampled_idxs_p1_begin,
                        tour_k_p1,
                        rng);
            const size_t p1_idx = *std::max_element(tour_sampled_idxs_p1_begin,
                                                    tour_sampled_idxs_p1.end(),
                                                    fitness_less_by_population_index);

            size_t p2_idx = p1_idx;
            if (population_size > 1)
            {
                std::sample(tour_idxs_begin,
                            tour_idxs_end - 1,
                            tour_sampled_idxs_p2_begin,
                            tour_k_p2,
                            rng);
                for (size_t j = 0; j < tour_k_p2; j++)
                {
                    size_t &sampled_idx = tour_sampled_idxs_p2[j];
                    if (sampled_idx >= p1_idx)
                        sampled_idx += 1;
                }
                p2_idx = *std::max_element(tour_sampled_idxs_p2_begin,
                                           tour_sampled_idxs_p2.end(),
                                           fitness_less_by_population_index);
            }

            const Individual &p1 = population[p1_idx];
            const Individual &p2 = population[p2_idx];

            Individual &child_1 = new_population[i];
            const bool has_second_child = (i + 1) < population_size;
            Individual &child_2 = has_second_child ? new_population[i + 1] : odd_aux_individual;

            if (prob_dist(rng) < GeneticPrivateConfig::CROSSOVER_PROBABILITY)
            {
                crossover_individual(p1,
                                     p2,
                                     child_1,
                                     child_2,
                                     has_second_child,
                                     ranges_12,
                                     ranges_21,
                                     rng);
            }
            else
            {
                child_1 = p1;
                if (has_second_child)
                    child_2 = p2;
            }

            mutate_individual(child_1, rng, ndist, prob_dist, mutation_probability_per_gene,
                              ranges_12, ranges_21);
            if (has_second_child)
                mutate_individual(child_2, rng, ndist, prob_dist, mutation_probability_per_gene,
                                  ranges_12, ranges_21);
        }

        std::swap(population, new_population);

        if (!calc_fitnesses(std::span<Individual>(population).subspan(ELITES),
                            fs,
                            effective_pad_12,
                            effective_pad_21,
                            pad_buffers))
        {
            QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
            return;
        }

        QMetaObject::invokeMethod(this, "set_generations_completed", Qt::QueuedConnection,
                                  Q_ARG(double, gen + 1));
    }

    const Individual &best = *std::min_element(population.begin(), population.end(), fitness_descending);

    const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);
    if (!wait_until_RT_read_idx_or_stop(curr_synapse_idx))
    {
        QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
        return;
    }

    const size_t new_synapse_idx = 1 - curr_synapse_idx;

    copy_individual_synapse_params_to_runtime(params_12[new_synapse_idx],
                                              best.params_12,
                                              use_i_fast_12,
                                              use_i_slow_12);
    copy_individual_synapse_params_to_runtime(params_21[new_synapse_idx],
                                              best.params_21,
                                              use_i_fast_21,
                                              use_i_slow_21);
    synapse_idx.store(new_synapse_idx, std::memory_order_release);

    QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
}
