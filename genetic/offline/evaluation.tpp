#include <algorithm>
#include <cmath>
using namespace kfr;

namespace EvaluationPrivateConfig
{
    static constexpr double PAD_LEN_FACTOR = 1.5;

    static constexpr int BUTTERWORTH_ORDER = 4;

    static constexpr double I_FAST_WEIGHT = 0.5;
    static constexpr double I_SLOW_WEIGHT = 0.5;

    static constexpr double EXPECTED_I_MARGIN_RANGE_FACTOR = 0.8;
}

namespace EvaluationPrivateConstants
{
    static constexpr double M_SLOW_INITIAL_VALUE = 0.0;
}

static double calculate_expected_i_dist_max(double expected_i_min,
                                            double expected_i_max)
{
    const double range = expected_i_max - expected_i_min;
    return (range + (range * BOPublicConfig::EXPECTED_I_MARGIN_FACTOR * 2.0)) *
           EvaluationPrivateConfig::EXPECTED_I_MARGIN_RANGE_FACTOR;
}

static double rescale_to_target(double value,
                                double src_min, double src_range,
                                double dst_min, double dst_range)
{
    const double norm = (value - src_min) / safe_divisor(src_range);
    return dst_min + (norm * dst_range);
}

static double rescale_to_target_no_offset(double value,
                                          double src_range,
                                          double dst_range)
{
    const double norm = value / safe_divisor(src_range);
    return norm * dst_range;
}

static double pearson_score(univector<double> &sig,
                            const univector<double> &ref_sig_centered,
                            double ref_sig_factor,
                            bool search_phase)
{
    sig -= mean(sig);
    const double sig_factor = std::sqrt(sum(sqr(sig)));
    const double r = sum(sig * ref_sig_centered) / safe_divisor(sig_factor * ref_sig_factor);
    const double normalized = (r + 1.0) / 2.0;
    return search_phase ? 1.0 - normalized : normalized;
}

static double range_score(double observed_min, double observed_max,
                          double expected_min, double expected_max,
                          double pts_dist_max)
{
    const double error = (std::abs(observed_min - expected_min) +
                          std::abs(observed_max - expected_max)) *
                         0.5;
    const double normalized_error = error / safe_divisor(pts_dist_max);
    return 1.0 - normalized_error;
}

ConstantEvaluationVals calc_constant_evaluation_vals(
    const univector_ref<const double> &v_pre_sig,
    double v_pre_min,
    double v_pre_max,
    double csv_step,
    double fc,
    bool use_i_fast,
    bool use_i_slow,
    double expected_i_min,
    double expected_i_max,
    bool search_phase)
{
    ConstantEvaluationVals result;

    const double fs = 1.0 / safe_divisor(csv_step);
    const size_t use_size = v_pre_sig.size();
    const size_t effective_pad = std::min(use_size - 1,
                                          static_cast<size_t>(EvaluationPrivateConfig::PAD_LEN_FACTOR * fs / safe_divisor(fc)));
    univector<double> padded(use_size + (2 * effective_pad));

    univector_ref<double> padded_seg = padded.slice(effective_pad, use_size);
    process(padded_seg, v_pre_sig);

    double *padded_ptr = padded.data();
    const double *v_pre_sig_ptr = v_pre_sig.data();
    const double left_edge_x2 = 2.0 * v_pre_sig_ptr[0];
    const double right_edge_x2 = 2.0 * v_pre_sig_ptr[use_size - 1];
    for (size_t i = 0; i < effective_pad; i++)
    {
        padded_ptr[effective_pad - 1 - i] = left_edge_x2 - v_pre_sig_ptr[i + 1];
        padded_ptr[use_size + effective_pad + i] = right_edge_x2 - v_pre_sig_ptr[use_size - 2 - i];
    }

    filtfilt(padded, to_sos<double>(iir_lowpass(
                              butterworth(EvaluationPrivateConfig::BUTTERWORTH_ORDER), fc, fs)));

    const bool use_both = use_i_fast && use_i_slow;
    double i_dist_max = 0.0;
    double v_pre_range = 0.0, expected_i_range = 0.0;
    if (use_both)
    {
        v_pre_range = v_pre_max - v_pre_min;
        expected_i_range = expected_i_max - expected_i_min;
    }
    else
    {
        i_dist_max = calculate_expected_i_dist_max(expected_i_min, expected_i_max);
    }

    if (use_i_fast)
    {
        univector<double> &ref_i_fast_sig_centered = result.ref_i_fast_sig_centered;
        ref_i_fast_sig_centered = v_pre_sig - padded_seg;

        double &ref_i_fast_min = result.ref_i_fast_min;
        double &ref_i_fast_max = result.ref_i_fast_max;
        double &i_fast_dist_max = result.i_fast_dist_max;
        if (use_both)
        {
            if (search_phase)
            {
                ref_i_fast_min = rescale_to_target_no_offset(-maxof(ref_i_fast_sig_centered), v_pre_range, expected_i_range);
                ref_i_fast_max = rescale_to_target_no_offset(-minof(ref_i_fast_sig_centered), v_pre_range, expected_i_range);
            }
            else
            {
                ref_i_fast_min = rescale_to_target_no_offset(minof(ref_i_fast_sig_centered), v_pre_range, expected_i_range);
                ref_i_fast_max = rescale_to_target_no_offset(maxof(ref_i_fast_sig_centered), v_pre_range, expected_i_range);
            }
            i_fast_dist_max = calculate_expected_i_dist_max(ref_i_fast_min, ref_i_fast_max);
        }
        else
        {
            ref_i_fast_min = expected_i_min;
            ref_i_fast_max = expected_i_max;
            i_fast_dist_max = i_dist_max;
        }

        ref_i_fast_sig_centered -= mean(ref_i_fast_sig_centered);
        result.ref_i_fast_sig_factor = std::sqrt(sum(sqr(ref_i_fast_sig_centered)));
    }

    if (use_i_slow)
    {
        univector<double> &ref_i_slow_sig_centered = result.ref_i_slow_sig_centered;
        ref_i_slow_sig_centered = padded_seg;

        double &ref_i_slow_min = result.ref_i_slow_min;
        double &ref_i_slow_max = result.ref_i_slow_max;
        double &i_slow_dist_max = result.i_slow_dist_max;
        if (use_both)
        {
            if (search_phase)
            {
                const double v_pre_min_to_use = -v_pre_max;
                ref_i_slow_min = rescale_to_target(-maxof(ref_i_slow_sig_centered), v_pre_min_to_use, v_pre_range,
                                                   expected_i_min, expected_i_range);
                ref_i_slow_max = rescale_to_target(-minof(ref_i_slow_sig_centered), v_pre_min_to_use, v_pre_range,
                                                   expected_i_min, expected_i_range);
            }
            else
            {
                ref_i_slow_min = rescale_to_target(minof(ref_i_slow_sig_centered), v_pre_min, v_pre_range,
                                                   expected_i_min, expected_i_range);
                ref_i_slow_max = rescale_to_target(maxof(ref_i_slow_sig_centered), v_pre_min, v_pre_range,
                                                   expected_i_min, expected_i_range);
            }
            i_slow_dist_max = calculate_expected_i_dist_max(ref_i_slow_min, ref_i_slow_max);
        }
        else
        {
            ref_i_slow_min = expected_i_min;
            ref_i_slow_max = expected_i_max;
            i_slow_dist_max = i_dist_max;
        }

        ref_i_slow_sig_centered -= mean(ref_i_slow_sig_centered);
        result.ref_i_slow_sig_factor = std::sqrt(sum(sqr(ref_i_slow_sig_centered)));
    }

    return result;
}

template <typename Integrator, typename NeuronType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
ChemicalSynapseEvaluation evaluate_candidate(
    const ChemicalSynapseParams &candidate,
    ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapse,
    NeuronType &model_neur,
    const ScaledSigResult &scaled_result,
    EvaluationISigBuffers &buffers,
    size_t v_pre_sig_start_idx,
    const ConstantEvaluationVals &constant_evaluation_vals,
    bool use_i_fast,
    bool use_i_slow,
    bool search_phase,
    double i_min,
    double i_max,
    ResetStateFuncType reset_state_neur,
    GetVFuncType get_v_neur)
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    const ChemicalSynapseParams &params = candidate;

    synapse.set(ChemicalSynapsisType::Esyn, params.e_syn);

    double *i_fast_sig_ptr = nullptr;
    if (use_i_fast)
    {
        i_fast_sig_ptr = buffers.i_fast_sig.data();

        synapse.set(ChemicalSynapsisType::gfast, params.g_fast);
        synapse.set(ChemicalSynapsisType::sfast, params.s_fast);
        synapse.set(ChemicalSynapsisType::Vfast, params.v_fast);
    }

    double *i_slow_sig_ptr = nullptr;
    if (use_i_slow)
    {
        i_slow_sig_ptr = buffers.i_slow_sig.data();

        synapse.set(ChemicalSynapsisType::gslow, params.g_slow);
        synapse.set(ChemicalSynapsisType::Vslow, params.v_slow);
        synapse.set(ChemicalSynapsisType::sslow, params.s_slow);
        synapse.set(ChemicalSynapsisType::k1, params.k1);
        synapse.set(ChemicalSynapsisType::k2, params.k2);
    }

    synapse.set(ChemicalSynapsisType::mslow, EvaluationPrivateConstants::M_SLOW_INITIAL_VALUE);

    reset_state_neur(model_neur);

    constexpr auto i_enum = ChemicalSynapsisType::i;
    constexpr auto i_fast_enum = ChemicalSynapsisType::ifast;
    constexpr auto i_slow_enum = ChemicalSynapsisType::islow;

    const size_t total_size = scaled_result.sig.size();
    const size_t points_factor = scaled_result.points_factor;
    const double dt = scaled_result.dt;
    const double *v_pre_sig_ptr = scaled_result.sig.data();
    const double *interpolated_points_ptr = scaled_result.interpolated_points.data();

    size_t interp_pts_counter = 0;
    size_t v_pre_sig_idx = 0;
    for (; v_pre_sig_idx < v_pre_sig_start_idx; v_pre_sig_idx++)
    {
        synapse.step(dt, v_pre_sig_ptr[v_pre_sig_idx], get_v_neur(model_neur));
        const double i_val = std::clamp(synapse.get(i_enum), i_min, i_max);
        model_neur.add_synaptic_input(-i_val);
        model_neur.step(dt);
        for (size_t k = 1; k < points_factor; k++)
        {
            synapse.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_neur(model_neur));
            const double i_interp_val = std::clamp(synapse.get(i_enum), i_min, i_max);
            model_neur.add_synaptic_input(-i_interp_val);
            model_neur.step(dt);
            interp_pts_counter++;
        }
    }

    size_t syn_sig_idx = 0;
    for (; v_pre_sig_idx < total_size - 1; v_pre_sig_idx++, syn_sig_idx++)
    {
        const double v_post = get_v_neur(model_neur);
        synapse.step(dt, v_pre_sig_ptr[v_pre_sig_idx], v_post);
        const double i_val = std::clamp(synapse.get(i_enum), i_min, i_max);
        model_neur.add_synaptic_input(-i_val);
        model_neur.step(dt);
        if (use_i_fast)
            i_fast_sig_ptr[syn_sig_idx] = synapse.get(i_fast_enum);
        if (use_i_slow)
            i_slow_sig_ptr[syn_sig_idx] = synapse.get(i_slow_enum);
        for (size_t k = 1; k < points_factor; k++)
        {
            synapse.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_neur(model_neur));
            const double i_interp_val = std::clamp(synapse.get(i_enum), i_min, i_max);
            model_neur.add_synaptic_input(-i_interp_val);
            model_neur.step(dt);
            interp_pts_counter++;
        }
    }

    const double v_post = get_v_neur(model_neur);
    synapse.step(dt, v_pre_sig_ptr[v_pre_sig_idx], v_post);
    const double i_val = std::clamp(synapse.get(i_enum), i_min, i_max);
    model_neur.add_synaptic_input(-i_val);
    model_neur.step(dt);
    if (use_i_fast)
        i_fast_sig_ptr[syn_sig_idx] = synapse.get(i_fast_enum);
    if (use_i_slow)
        i_slow_sig_ptr[syn_sig_idx] = synapse.get(i_slow_enum);

    constexpr double I_FAST_WEIGHT = EvaluationPrivateConfig::I_FAST_WEIGHT;
    constexpr double I_SLOW_WEIGHT = EvaluationPrivateConfig::I_SLOW_WEIGHT;

    double i_range_score_accum = 0.0;
    double i_shape_score_accum = 0.0;
    double total_weight = 0.0;

    if (use_i_fast)
    {
        univector<double> &i_fast_sig = buffers.i_fast_sig;
        const double i_fast_range_score = range_score(
            minof(i_fast_sig),
            maxof(i_fast_sig),
            constant_evaluation_vals.ref_i_fast_min,
            constant_evaluation_vals.ref_i_fast_max,
            constant_evaluation_vals.i_fast_dist_max);

        const double i_fast_shape_score = pearson_score(
            i_fast_sig,
            constant_evaluation_vals.ref_i_fast_sig_centered,
            constant_evaluation_vals.ref_i_fast_sig_factor,
            search_phase);

        i_range_score_accum += I_FAST_WEIGHT * i_fast_range_score;
        i_shape_score_accum += I_FAST_WEIGHT * i_fast_shape_score;
        total_weight += I_FAST_WEIGHT;
    }

    if (use_i_slow)
    {
        univector<double> &i_slow_sig = buffers.i_slow_sig;
        const double i_slow_range_score = range_score(
            minof(i_slow_sig),
            maxof(i_slow_sig),
            constant_evaluation_vals.ref_i_slow_min,
            constant_evaluation_vals.ref_i_slow_max,
            constant_evaluation_vals.i_slow_dist_max);

        const double i_slow_shape_score = pearson_score(
            i_slow_sig,
            constant_evaluation_vals.ref_i_slow_sig_centered,
            constant_evaluation_vals.ref_i_slow_sig_factor,
            search_phase);

        i_range_score_accum += I_SLOW_WEIGHT * i_slow_range_score;
        i_shape_score_accum += I_SLOW_WEIGHT * i_slow_shape_score;
        total_weight += I_SLOW_WEIGHT;
    }

    return ChemicalSynapseEvaluation(i_range_score_accum / total_weight,
                                     i_shape_score_accum / total_weight);
}
