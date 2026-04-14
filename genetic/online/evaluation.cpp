
#include <chrono>
#include <cmath>
#include <kfr/all.hpp>
#include "bidirectional_chemical_synapse_BO.h"
#include <concepts>
using namespace kfr;

namespace EvaluationPrivateConfig
{
    // Butterworth low-pass filter order used by filtfilt (zero-phase) separation.
    static constexpr int BUTTERWORTH_ORDER = 4;

    // Relative weights for the fast vs slow components when both are enabled.
    static constexpr double I_FAST_WEIGHT = 0.5;
    static constexpr double I_SLOW_WEIGHT = 0.5;
}

// Linear mapping from a source interval [src_min, src_min+src_range]
// to a destination interval [dst_min, dst_min+dst_range].
static double rescale_to_target(double value,
                                double src_min, double src_range,
                                double dst_min, double dst_range)
{
    const double norm = (value - src_min) / safe_divisor(src_range);
    return dst_min + (norm * dst_range);
}

template <typename T>
    requires std::same_as<T, univector<double>> || std::same_as<T, univector_ref<double>>
static double pearson_score(univector<double> &sig,
                            const T &ref_sig_centered,
                            double ref_sig_factor,
                            unsigned int search_phase)
{
    // Pearson correlation based shape similarity.
    // Output is normalized to [0,1]. `search_phase` optionally flips the target.
    sig -= mean(sig);
    const double sig_factor = std::sqrt(sum(sqr(sig)));
    const double r = sum(sig * ref_sig_centered) / safe_divisor(sig_factor * ref_sig_factor);
    const double normalized = (r + 1.0) / 2.0;
    return search_phase ? 1.0 - normalized : normalized;
}

// Scores how close the observed (min,max) is to the expected (min,max).
static double range_score(double observed_min, double observed_max,
                          double expected_min, double expected_max,
                          double max_pts_dist)
{
    const double error = (std::abs(observed_min - expected_min) +
                          std::abs(observed_max - expected_max)) *
                         0.5;
    const double normalized_error = error / safe_divisor(max_pts_dist);
    return 1.0 - normalized_error;
}

static ChemicalSynapseEvaluation evaluate_sigs_one_direction(
    univector<double> &vpre_sig,
    univector<double> &i_fast_sig,
    univector<double> &i_slow_sig,
    size_t effective_pad,
    double fs,
    double fc,
    unsigned int use_i_fast,
    unsigned int use_i_slow,
    unsigned int search_phase,
    double expected_i_min,
    double expected_i_max,
    double max_i_dist,
    univector<double> &padded_buff)
{
    // Pipeline:
    // 1) Copy vpre into a padded buffer and mirror-pad edges (reduces filter transients)
    // 2) Zero-phase low-pass filter (filtfilt) -> reference slow component
    // 3) Reference fast component = (raw vpre) - (low-pass vpre)
    // 4) Compute range + shape scores for enabled currents
    const size_t use_size = vpre_sig.size();

    univector_ref<double> padded_seg = padded_buff.slice(effective_pad, use_size);
    process(padded_seg, vpre_sig);

    double *padded_ptr = padded_buff.data();
    const double *vpre_sig_ptr = vpre_sig.data();
    const double left_edge_x2 = 2.0 * vpre_sig_ptr[0];
    const double right_edge_x2 = 2.0 * vpre_sig_ptr[use_size - 1];
    for (size_t i = 0; i < effective_pad; i++)
    {
        // Mirror padding around the endpoints: x[-i] = 2*x0 - x[i]
        // This keeps the first derivative closer to continuous at the edges.
        padded_ptr[effective_pad - 1 - i] = left_edge_x2 - vpre_sig_ptr[i + 1];
        padded_ptr[use_size + effective_pad + i] = right_edge_x2 - vpre_sig_ptr[use_size - 2 - i];
    }

    // Zero-phase low-pass (slow component reference).
    filtfilt(padded_buff, to_sos<double>(iir_lowpass(
                              butterworth(EvaluationPrivateConfig::BUTTERWORTH_ORDER), fc, fs)));

    const bool use_both = use_i_fast && use_i_slow;

    constexpr double I_FAST_WEIGHT = EvaluationPrivateConfig::I_FAST_WEIGHT;
    constexpr double I_SLOW_WEIGHT = EvaluationPrivateConfig::I_SLOW_WEIGHT;

    double i_range_score_accum = 0.0;
    double i_shape_score_accum = 0.0;
    double total_weight = 0.0;

    double vpre_min = 0.0, vpre_range = 0.0, expected_i_range = 0.0;
    if (use_both)
    {
        // When both components are enabled, we derive per-component expected ranges
        // by rescaling reference signals into the user-provided expected_i interval.
        vpre_min = minof(vpre_sig);
        const double vpre_max = maxof(vpre_sig);
        vpre_range = vpre_max - vpre_min;
        expected_i_range = expected_i_max - expected_i_min;
    }

    if (use_i_fast)
    {
        // NOTE: vpre_sig is reused as a temporary buffer to avoid allocations.
        // ref_i_fast_sig_aux := (raw_vpre - lowpass_vpre) -> high-frequency reference.
        univector<double> &ref_i_fast_sig_aux = vpre_sig;
        ref_i_fast_sig_aux -= padded_seg;

        double ref_i_fast_min, ref_i_fast_max, max_i_fast_dist;
        if (use_both)
        {
            // Rescale reference fast min/max into expected current bounds.
            ref_i_fast_min = rescale_to_target(minof(ref_i_fast_sig_aux), vpre_min, vpre_range,
                                               expected_i_min, expected_i_range);
            ref_i_fast_max = rescale_to_target(maxof(ref_i_fast_sig_aux), vpre_min, vpre_range,
                                               expected_i_min, expected_i_range);
            max_i_fast_dist = calculate_expected_i_max_dist(ref_i_fast_min, ref_i_fast_max);
        }
        else
        {
            ref_i_fast_min = expected_i_min;
            ref_i_fast_max = expected_i_max;
            max_i_fast_dist = max_i_dist;
        }

        ref_i_fast_sig_aux -= mean(ref_i_fast_sig_aux);
        const double ref_i_fast_factor = std::sqrt(sum(sqr(ref_i_fast_sig_aux)));

        const double i_fast_range_score = range_score(
            minof(i_fast_sig),
            maxof(i_fast_sig),
            ref_i_fast_min,
            ref_i_fast_max,
            max_i_fast_dist);

        const double i_fast_shape_score =
            pearson_score(i_fast_sig,
                          ref_i_fast_sig_aux,
                          ref_i_fast_factor,
                          search_phase);

        i_range_score_accum += I_FAST_WEIGHT * i_fast_range_score;
        i_shape_score_accum += I_FAST_WEIGHT * i_fast_shape_score;
        total_weight += I_FAST_WEIGHT;
    }

    if (use_i_slow)
    {
        // Slow reference is the low-pass filtered vpre segment.
        univector_ref<double> &ref_i_slow_sig_aux = padded_seg;

        double ref_i_slow_min, ref_i_slow_max, max_i_slow_dist;
        if (use_both)
        {
            // Rescale reference slow min/max into expected current bounds.
            ref_i_slow_min = rescale_to_target(minof(ref_i_slow_sig_aux), vpre_min, vpre_range,
                                               expected_i_min, expected_i_range);
            ref_i_slow_max = rescale_to_target(maxof(ref_i_slow_sig_aux), vpre_min, vpre_range,
                                               expected_i_min, expected_i_range);
            max_i_slow_dist = calculate_expected_i_max_dist(ref_i_slow_min, ref_i_slow_max);
        }
        else
        {
            ref_i_slow_min = expected_i_min;
            ref_i_slow_max = expected_i_max;
            max_i_slow_dist = max_i_dist;
        }

        ref_i_slow_sig_aux -= mean(ref_i_slow_sig_aux);
        const double ref_i_slow_factor = std::sqrt(sum(sqr(ref_i_slow_sig_aux)));

        const double i_slow_range_score = range_score(
            minof(i_slow_sig),
            maxof(i_slow_sig),
            ref_i_slow_min,
            ref_i_slow_max,
            max_i_slow_dist);

        const double i_slow_shape_score =
            pearson_score(i_slow_sig,
                          ref_i_slow_sig_aux,
                          ref_i_slow_factor,
                          search_phase);

        i_range_score_accum += I_SLOW_WEIGHT * i_slow_range_score;
        i_shape_score_accum += I_SLOW_WEIGHT * i_slow_shape_score;
        total_weight += I_SLOW_WEIGHT;
    }

    ChemicalSynapseEvaluation result;
    result.i_range_score = i_range_score_accum / total_weight;
    result.i_shape_score = i_shape_score_accum / total_weight;
    return result;
}

ChemicalSynapseEvaluation BidirectionalChemicalSynapseBO::evaluate_candidate(
    const Candidate &candidate,
    double fs,
    size_t effective_pad_12,
    size_t effective_pad_21,
    EvaluationPadBuffers &pad_buffers,
    double max_i_dist_12,
    double max_i_dist_21,
    size_t &curr_synapse_idx)
{
    // Evaluation is called from the BO thread. It must cooperate with the RT loop.
    if (stop_BO.load(std::memory_order_relaxed))
        throw StopEvaluation();

    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

    if (!use_syn_12 && !use_syn_21)
    {
        // Nothing to optimize: still advance the evaluation counter for UI consistency.
        QMetaObject::invokeMethod(this, "set_evaluations_completed", Qt::QueuedConnection,
                                  Q_ARG(double, evaluations_completed + 1));
        return ChemicalSynapseEvaluation{0.0, 0.0};
    }

    const std::chrono::duration<double, std::milli> active_wait_duration(BOPublicConfig::ACTIVE_WAIT_MS);

    // Wait until RT has read the current slot so we can safely swap to the other.
    if (!wait_until_RT_read_idx_or_stop(curr_synapse_idx))
        throw StopEvaluation();

    const size_t new_synapse_idx = 1 - curr_synapse_idx;

    // Write candidate params into inactive slot (NRT-only write).
    copy_selected_synapse_params(params_12[new_synapse_idx],
                                 candidate.params_12,
                                 use_i_fast_12,
                                 use_i_slow_12);

    copy_selected_synapse_params(params_21[new_synapse_idx],
                                 candidate.params_21,
                                 use_i_fast_21,
                                 use_i_slow_21);

    // Publish new params to RT (release pairs with RT acquire-load in execute()).
    synapse_idx.store(new_synapse_idx, std::memory_order_release);
    curr_synapse_idx = new_synapse_idx;

    if (stabilization_time > 0.0)
    {
        std::this_thread::sleep_for(std::chrono::duration<double, std::milli>(stabilization_time));
    }

    storing_idx = 0;

    // Ask RT thread to start capturing signals into pre-allocated buffers.
    RT_storing.store(true, std::memory_order_release);

    // Active-wait until RT indicates capture finished (or stop requested).
    while (RT_storing.load(std::memory_order_acquire))
    {
        if (stop_BO.load(std::memory_order_relaxed))
        {
            RT_storing.store(false, std::memory_order_relaxed);
            throw StopEvaluation();
        }
        std::this_thread::sleep_for(active_wait_duration);
    }

    double i_range_score_accum = 0.0;
    double i_shape_score_accum = 0.0;
    double num_directions = 0.0;

    if (use_syn_12)
    {
        // Score direction 1->2.
        ChemicalSynapseEvaluation score_12 = evaluate_sigs_one_direction(
            v_sig_1,
            i_fast_sig_12,
            i_slow_sig_12,
            effective_pad_12,
            fs,
            fc_1,
            use_i_fast_12,
            use_i_slow_12,
            search_phase,
            expected_i_min_12, expected_i_max_12,
            max_i_dist_12,
            pad_buffers.padded_buff_12);
        i_range_score_accum += score_12.i_range_score;
        i_shape_score_accum += score_12.i_shape_score;
        num_directions++;
    }

    if (use_syn_21)
    {
        // Score direction 2->1.
        ChemicalSynapseEvaluation score_21 = evaluate_sigs_one_direction(
            v_sig_2,
            i_fast_sig_21,
            i_slow_sig_21,
            effective_pad_21,
            fs,
            fc_2,
            use_i_fast_21,
            use_i_slow_21,
            search_phase,
            expected_i_min_21, expected_i_max_21,
            max_i_dist_21,
            pad_buffers.padded_buff_21);
        i_range_score_accum += score_21.i_range_score;
        i_shape_score_accum += score_21.i_shape_score;
        num_directions++;
    }

    const double i_range_score = i_range_score_accum / num_directions;
    const double i_shape_score = i_shape_score_accum / num_directions;
    ChemicalSynapseEvaluation result;
    // Avoid NaN/Inf propagating into BO.
    result.i_range_score = std::isfinite(i_range_score) ? i_range_score : 0.0;
    result.i_shape_score = std::isfinite(i_shape_score) ? i_shape_score : 0.0;

    QMetaObject::invokeMethod(this, "set_evaluations_completed", Qt::QueuedConnection,
                              Q_ARG(double, evaluations_completed + 1));

    return result;
}
