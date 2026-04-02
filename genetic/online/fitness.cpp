
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <kfr/all.hpp>
#include "bidirectional_chemical_synapse_genetic.h"
#include "utils.hpp"
#include <concepts>
using namespace kfr;

namespace FitnessPrivateConfig
{
    static constexpr double SYN_RANGE_WEIGHT = 0.5;
    static constexpr double VPRE_SYN_COMP_WEIGHT = 0.5;

    static constexpr double IFAST_WEIGHT = 0.5;
    static constexpr double ISLOW_WEIGHT = 0.5;

    static constexpr int BUTTERWORTH_ORDER = 4;
}

static inline double rescale_to_target(double value,
                                       double src_min, double src_range,
                                       double dst_min, double dst_range)
{
    const double norm = (value - src_min) / src_range;
    return dst_min + (norm * dst_range);
}

template <typename T>
    requires std::same_as<T, univector<double>> || std::same_as<T, univector_ref<double>>
static inline double pearson_score(univector<double> &sig,
                                   const T &ref_sig_centered,
                                   double ref_sig_factor,
                                   unsigned int search_phase)
{
    sig -= mean(sig);
    const double sig_factor = std::sqrt(sum(sqr(sig)));
    const double r = sum(sig * ref_sig_centered) / (sig_factor * ref_sig_factor);
    const double normalized = (r + 1.0) / 2.0;
    return search_phase ? 1.0 - normalized : normalized;
}

static inline double calc_range_score_component(double observed_min, double observed_max,
                                                double expected_min, double expected_max,
                                                double max_diff)
{
    return 1.0 - (((std::abs(observed_min - expected_min) +
                    std::abs(observed_max - expected_max)) *
                   0.5) /
                  max_diff);
}

static double calc_fitness_from_sigs_one_direction(
    univector<double> &vpre_sig,
    univector<double> &i_fast_sig,
    univector<double> &i_slow_sig,
    size_t effective_pad,
    double fs,
    unsigned int use_i_fast,
    unsigned int use_i_slow,
    unsigned int search_phase,
    double expected_i_min,
    double expected_i_max,
    univector<double> &padded_buff)
{
    const size_t use_size = vpre_sig.size();

    univector_ref<double> padded_seg = padded_buff.slice(effective_pad, use_size);
    process(padded_seg, vpre_sig);

    double *padded_ptr = padded_buff.data();
    const double *vpre_sig_ptr = vpre_sig.data();
    const double left_edge_x2 = 2.0 * vpre_sig_ptr[0];
    const double right_edge_x2 = 2.0 * vpre_sig_ptr[use_size - 1];
    for (size_t i = 0; i < effective_pad; i++)
    {
        padded_ptr[effective_pad - 1 - i] = left_edge_x2 - vpre_sig_ptr[i + 1];
        padded_ptr[use_size + effective_pad + i] = right_edge_x2 - vpre_sig_ptr[use_size - 2 - i];
    }

    filtfilt(padded_buff, to_sos<double>(iir_lowpass(
                              butterworth(FitnessPrivateConfig::BUTTERWORTH_ORDER), FitnessPublicConfig::FILTER_FC, fs)));

    const bool use_both = use_i_fast && use_i_slow;

    constexpr double IFAST_WEIGHT = FitnessPrivateConfig::IFAST_WEIGHT;
    constexpr double ISLOW_WEIGHT = FitnessPrivateConfig::ISLOW_WEIGHT;

    double vpre_syn_comp_accum = 0.0;
    double i_range_accum = 0.0;
    double total_weight = 0.0;

    double vpre_min = 0.0, vpre_range = 0.0, expected_i_range = 0.0;
    if (use_both)
    {
        vpre_min = minof(vpre_sig);
        const double vpre_max = maxof(vpre_sig);
        vpre_range = vpre_max - vpre_min;
        expected_i_range = expected_i_max - expected_i_min;
    }

    if (use_i_fast)
    {
        univector<double> &ref_i_fast_sig_aux = vpre_sig;
        ref_i_fast_sig_aux -= padded_seg;

        double ref_i_fast_min, ref_i_fast_max;
        if (use_both)
        {
            ref_i_fast_min = rescale_to_target(minof(ref_i_fast_sig_aux), vpre_min, vpre_range,
                                               expected_i_min, expected_i_range);
            ref_i_fast_max = rescale_to_target(maxof(ref_i_fast_sig_aux), vpre_min, vpre_range,
                                               expected_i_min, expected_i_range);
        }
        else
        {
            ref_i_fast_min = expected_i_min;
            ref_i_fast_max = expected_i_max;
        }

        ref_i_fast_sig_aux -= mean(ref_i_fast_sig_aux);
        const double ref_i_fast_factor = std::sqrt(sum(sqr(ref_i_fast_sig_aux)));

        const double comp_ifast = calc_range_score_component(
            minof(i_fast_sig),
            maxof(i_fast_sig),
            ref_i_fast_min,
            ref_i_fast_max,
            ref_i_fast_max - ref_i_fast_min);
        i_range_accum += IFAST_WEIGHT * comp_ifast;

        const double vpost_ifast_comp_score =
            pearson_score(i_fast_sig,
                          ref_i_fast_sig_aux,
                          ref_i_fast_factor,
                          search_phase);
        vpre_syn_comp_accum += IFAST_WEIGHT * vpost_ifast_comp_score;

        total_weight += IFAST_WEIGHT;
    }

    if (use_i_slow)
    {
        univector_ref<double> &ref_i_slow_sig_aux = padded_seg;

        double ref_i_slow_min, ref_i_slow_max;
        if (use_both)
        {
            ref_i_slow_min = rescale_to_target(minof(ref_i_slow_sig_aux), vpre_min, vpre_range,
                                               expected_i_min, expected_i_range);
            ref_i_slow_max = rescale_to_target(maxof(ref_i_slow_sig_aux), vpre_min, vpre_range,
                                               expected_i_min, expected_i_range);
        }
        else
        {
            ref_i_slow_min = expected_i_min;
            ref_i_slow_max = expected_i_max;
        }

        ref_i_slow_sig_aux -= mean(ref_i_slow_sig_aux);
        const double ref_i_slow_factor = std::sqrt(sum(sqr(ref_i_slow_sig_aux)));

        const double comp_islow = calc_range_score_component(
            minof(i_slow_sig),
            maxof(i_slow_sig),
            ref_i_slow_min,
            ref_i_slow_max,
            ref_i_slow_max - ref_i_slow_min);
        i_range_accum += ISLOW_WEIGHT * comp_islow;

        const double vpost_islow_comp_score =
            pearson_score(i_slow_sig,
                          ref_i_slow_sig_aux,
                          ref_i_slow_factor,
                          search_phase);
        vpre_syn_comp_accum += ISLOW_WEIGHT * vpost_islow_comp_score;

        total_weight += ISLOW_WEIGHT;
    }

    if (total_weight == 0.0)
        return 0.0;

    return ((FitnessPrivateConfig::SYN_RANGE_WEIGHT * i_range_accum) +
            (FitnessPrivateConfig::VPRE_SYN_COMP_WEIGHT * vpre_syn_comp_accum)) /
           total_weight;
}

double BidirectionalChemicalSynapseGenetic::calc_fitness_from_sigs(double fs,
                                                                   size_t effective_pad,
                                                                   univector<double> &padded_buff)
{
    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

    double score_12 = 0.0;
    double score_21 = 0.0;
    double num_directions = 0.0;

    if (use_syn_12)
    {
        score_12 = calc_fitness_from_sigs_one_direction(
            v_sig_1,
            i_fast_sig_12,
            i_slow_sig_12,
            effective_pad,
            fs,
            use_i_fast_12,
            use_i_slow_12,
            search_phase,
            expected_i_min_12, expected_i_max_12,
            padded_buff);

        num_directions++;
    }

    if (use_syn_21)
    {
        score_21 = calc_fitness_from_sigs_one_direction(
            v_sig_2,
            i_fast_sig_21,
            i_slow_sig_21,
            effective_pad,
            fs,
            use_i_fast_21,
            use_i_slow_21,
            search_phase,
            expected_i_min_21, expected_i_max_21,
            padded_buff);

        num_directions++;
    }

    if (num_directions == 0.0)
        return 0.0;

    const double final_score = (score_12 + score_21) / num_directions;

    if (!(std::isfinite(final_score)) || final_score < 0.0)
        return 0.0;

    return final_score;
}

bool BidirectionalChemicalSynapseGenetic::calc_fitnesses(std::span<Individual> individuals,
                                                         double fs,
                                                         size_t effective_pad,
                                                         univector<double> &padded_buff)
{
    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

    if (!use_syn_12 && !use_syn_21)
    {
        for (Individual &individual : individuals)
        {
            individual.fitness = 0.0;
        }
        return true;
    }

    const std::chrono::duration<double, std::milli> stabilization_duration(stabilization_time);
    const std::chrono::duration<double, std::milli> active_wait_duration(GeneticPublicConfig::ACTIVE_WAIT_MS);
    size_t curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);

    for (Individual &individual : individuals)
    {
        if (stop_genetic.load(std::memory_order_relaxed))
            return false;

        if (!wait_until_RT_read_idx_or_stop(curr_synapse_idx))
            return false;

        const size_t new_synapse_idx = 1 - curr_synapse_idx;

        if (use_syn_12)
        {
            apply_variation_params(params_12[new_synapse_idx],
                                   individual.variation_params_12,
                                   use_i_fast_12,
                                   use_i_slow_12);
        }

        if (use_syn_21)
        {
            apply_variation_params(params_21[new_synapse_idx],
                                   individual.variation_params_21,
                                   use_i_fast_21,
                                   use_i_slow_21);
        }

        synapse_idx.store(new_synapse_idx, std::memory_order_release);
        curr_synapse_idx = new_synapse_idx;

        if (stabilization_time > 0.0)
        {
            std::this_thread::sleep_for(stabilization_duration);
        }

        storing_idx = 0;

        RT_storing.store(true, std::memory_order_release);

        while (RT_storing.load(std::memory_order_acquire))
        {
            if (stop_genetic.load(std::memory_order_relaxed))
            {
                RT_storing.store(false, std::memory_order_relaxed);
                return false;
            }
            std::this_thread::sleep_for(active_wait_duration);
        }

        individual.fitness = calc_fitness_from_sigs(fs, effective_pad, padded_buff);

        QMetaObject::invokeMethod(this, "set_individuals_completed", Qt::QueuedConnection,
                                  Q_ARG(double, individuals_completed + 1));
    }

    return true;
}
