
#include <algorithm>
#include <cmath>
#include <iostream>
#include <kfr/all.hpp>

namespace FitnessConfig
{
    static constexpr double V_COMP_WEIGHT = 0.5;
    static constexpr double VPRE_I_COMP_WEIGHT = 0.5;
    static constexpr double VPOST_IFAST_COMP_WEIGHT = 0.5;
    static constexpr double VPOST_ISLOW_COMP_WEIGHT = 0.5;

    static constexpr double AVG_SMOOTH_POINTS_BURST_DIVISOR = 100;

    static constexpr size_t FILTER_PAD_LEN = 1000;
    static constexpr int BUTTERWORTH_ORDER = 4;
    static constexpr double FILTER_FC_BURST_FRACTION = 0.5;
}

namespace FitnessConstants
{
    static constexpr double M_SLOW_INITIAL_VALUE = 0.0;

    static constexpr size_t MIN_AVG_SMOOTH_POINTS = 1;
}

ConstantSignalFitnessVals calc_const_signal_fitness_vals(const kfr::univector<double> &signal, double min_val, double max_val, bool search_phase, size_t avg_smooth_points, size_t start_index, double fs, double fc, SignalBuffers &buffers, bool use_ifast, bool use_islow)
{
    ConstantSignalFitnessVals result;

    const size_t use_size = signal.size() - start_index;

    kfr::univector<double> &smoothed_signal_to_fit = result.smoothed_signal_to_fit;
    smoothed_signal_to_fit.resize(use_size);
    kfr::univector<double> &normalized_signal_to_fit = result.normalized_signal_to_fit;
    normalized_signal_to_fit.resize(use_size);

    double running_sum = 0.0;
    double smoothed_min_val = GeneralConstants::DOUBLE_MAX;
    double smoothed_max_val = GeneralConstants::DOUBLE_MIN;

    const kfr::univector<double> signal_seg = signal.slice(start_index, start_index + use_size);
    const kfr::univector<double> prefix_seg = signal.slice(start_index - avg_smooth_points, start_index);

    running_sum = kfr::sum(prefix_seg);

    for (size_t sig_i = start_index, res_i = 0; res_i < use_size; res_i++, sig_i++)
    {
        const double sig_val = signal[sig_i];
        running_sum += sig_val;
        running_sum -= signal[sig_i - avg_smooth_points];
        const double smoothed_val = running_sum / avg_smooth_points;
        smoothed_signal_to_fit[res_i] = smoothed_val;
    }

    smoothed_min_val = kfr::minof(smoothed_signal_to_fit);
    smoothed_max_val = kfr::maxof(smoothed_signal_to_fit);

    normalized_signal_to_fit = (signal_seg - min_val) / (max_val - min_val);

    result.max_v_comp_distance = use_size * (smoothed_max_val - smoothed_min_val);

    calc_ifast_islow_ref_signals(signal, start_index, use_size, result, fs, fc, buffers, use_ifast, use_islow);

    if (!search_phase)
    {
        const double smoothed_flip_offset = smoothed_min_val + smoothed_max_val;
        smoothed_signal_to_fit = smoothed_flip_offset - smoothed_signal_to_fit;

        normalized_signal_to_fit = 1.0 - normalized_signal_to_fit;

        if (use_islow)
        {
            kfr::univector<double> &norm_signal_to_fit_islow = result.norm_signal_to_fit_islow;
            norm_signal_to_fit_islow = 1.0 - norm_signal_to_fit_islow;
        }
        if (use_ifast)
        {
            kfr::univector<double> &norm_signal_to_fit_ifast = result.norm_signal_to_fit_ifast;
            norm_signal_to_fit_ifast = 1.0 - norm_signal_to_fit_ifast;
        }
    }

    return result;
}

static void normalize_inplace(kfr::univector<double> &signal)
{
    const double sig_min = kfr::minof(signal);
    const double sig_max = kfr::maxof(signal);

    signal = (signal - sig_min) / (sig_max - sig_min);
}

static double compute_norm_component_score(const kfr::univector<double> &signal,
                                           const kfr::univector<double> &ref_norm)
{
    const size_t use_size = signal.size();

    const double min_val = kfr::minof(signal);
    const double max_val = kfr::maxof(signal);

    kfr::univector<double> norm = (signal - min_val) / (max_val - min_val);

    const double comp_dist = kfr::sum(kfr::abs(ref_norm - norm));

    return 1.0 - (comp_dist / use_size);
}

static void calc_ifast_islow_ref_signals(const kfr::univector<double> &signal,
                                         size_t start_index, size_t use_size,
                                         ConstantSignalFitnessVals &result,
                                         double fs, double fc,
                                         SignalBuffers &buffers,
                                         bool use_ifast, bool use_islow)
{
    kfr::univector<double> &padded = buffers.kfr_padded;

    constexpr size_t FILTER_PAD_LEN = FitnessConfig::FILTER_PAD_LEN;

    const auto padded_seg_begin = padded.begin() + FILTER_PAD_LEN;
    const auto signal_seg_begin = signal.cbegin() + start_index;
    const auto signal_seg_end = signal_seg_begin + use_size;

    std::copy(signal_seg_begin, signal_seg_end, padded_seg_begin);
    for (size_t i = 0; i < FILTER_PAD_LEN; i++)
    {
        padded[FILTER_PAD_LEN - 1 - i] = signal_seg_begin[i + 1];
        padded[use_size + FILTER_PAD_LEN + i] = signal_seg_begin[use_size - 2 - i];
    }

    kfr::filtfilt(padded, kfr::to_sos<double>(kfr::iir_lowpass(
                              kfr::butterworth(FitnessConfig::BUTTERWORTH_ORDER), fc, fs)));

    if (use_islow)
    {
        kfr::univector<double> &norm_signal_to_fit_islow = result.norm_signal_to_fit_islow;
        norm_signal_to_fit_islow.resize(use_size);
        std::copy(padded_seg_begin, padded_seg_begin + use_size, norm_signal_to_fit_islow.begin());
        normalize_inplace(norm_signal_to_fit_islow);
    }

    if (use_ifast)
    {
        kfr::univector<double> &norm_signal_to_fit_ifast = result.norm_signal_to_fit_ifast;
        norm_signal_to_fit_ifast.resize(use_size);
        std::copy(signal_seg_begin, signal_seg_end, norm_signal_to_fit_ifast.begin());
        norm_signal_to_fit_ifast -= padded.slice(FILTER_PAD_LEN, FILTER_PAD_LEN + use_size);
        normalize_inplace(norm_signal_to_fit_ifast);
    }
}

static double fitness_from_signals(const ConstantSignalFitnessVals &living_const_signal_fitness_vals, bool search_phase, size_t avg_smooth_points, bool use_ifast, bool use_islow, const SignalBuffers &buffers)
{
    const kfr::univector<double> &model_signal = buffers.model_signal;
    const kfr::univector<double> &synapsis_signal = buffers.synapsis_signal;

    double running_sum = 0.0;
    for (size_t i = 0; i < avg_smooth_points; i++)
        running_sum += model_signal[i];

    const size_t use_size = synapsis_signal.size();

    double v_comp_dist = 0.0;
    const kfr::univector<double> &living_smoothed_signal_to_fit = living_const_signal_fitness_vals.smoothed_signal_to_fit;
    for (size_t i = 0; i < use_size; i++)
    {
        running_sum += model_signal[i + avg_smooth_points];
        running_sum -= model_signal[i];
        const double smoothed_val = running_sum / avg_smooth_points;
        v_comp_dist += std::abs(living_smoothed_signal_to_fit[i] - smoothed_val);
    }
    const double v_comp_score = 1.0 - (v_comp_dist / living_const_signal_fitness_vals.max_v_comp_distance);

    double weighted_sum = FitnessConfig::V_COMP_WEIGHT * v_comp_score;
    double total_weight = FitnessConfig::V_COMP_WEIGHT;

    if (use_ifast && use_islow)
    {
        const kfr::univector<double> &living_norm_signal_to_fit =
            living_const_signal_fitness_vals.normalized_signal_to_fit;
        const double vpre_i_comp_score =
            compute_norm_component_score(buffers.synapsis_signal,
                                         living_norm_signal_to_fit);

        weighted_sum += FitnessConfig::VPRE_I_COMP_WEIGHT * vpre_i_comp_score;
        total_weight += FitnessConfig::VPRE_I_COMP_WEIGHT;
    }

    if (use_ifast)
    {
        const double vpost_ifast_comp_score =
            compute_norm_component_score(buffers.ifast_signal,
                                         living_const_signal_fitness_vals.norm_signal_to_fit_ifast);
        weighted_sum += FitnessConfig::VPOST_IFAST_COMP_WEIGHT * vpost_ifast_comp_score;
        total_weight += FitnessConfig::VPOST_IFAST_COMP_WEIGHT;
    }

    if (use_islow)
    {
        const double vpost_islow_comp_score =
            compute_norm_component_score(buffers.islow_signal,
                                         living_const_signal_fitness_vals.norm_signal_to_fit_islow);
        weighted_sum += FitnessConfig::VPOST_ISLOW_COMP_WEIGHT * vpost_islow_comp_score;
        total_weight += FitnessConfig::VPOST_ISLOW_COMP_WEIGHT;
    }

    const double final_score = weighted_sum / total_weight;

    if (std::isnan(final_score) || final_score < 0.0)
        return 0.0;

    return final_score;
}

template <typename Integrator, typename NeuronType, size_t N, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    std::array<Individual, N> &individuals,
                    const ScaledSignalResult &scaled_result,
                    const ConstantSignalFitnessVals &living_const_signal_fitness_vals,
                    bool search_phase,
                    SignalBuffers &buffers,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur,
                    size_t ind_start_index,
                    size_t signal_start_index,
                    size_t avg_smooth_points,
                    bool use_ifast,
                    bool use_islow)
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    constexpr auto i_enum = ChemicalSynapsisType::i;
    constexpr auto ifast_enum = ChemicalSynapsisType::ifast;
    constexpr auto islow_enum = ChemicalSynapsisType::islow;

    kfr::univector<double> &model_signal_buffer = buffers.model_signal;
    kfr::univector<double> &synapsis_signal_buffer = buffers.synapsis_signal;
    kfr::univector<double> &ifast_signal_buffer = buffers.ifast_signal;
    kfr::univector<double> &islow_signal_buffer = buffers.islow_signal;

    const size_t total_size = scaled_result.signal.size();
    const size_t points_factor = scaled_result.points_factor;
    const double dt = scaled_result.dt;
    const double *signal_data = scaled_result.signal.data();
    const double *interpolated_signal_data = scaled_result.interpolated_points.data();

    const size_t sig_start_minus_smoothed_avg_pts = signal_start_index - avg_smooth_points;

    for (size_t i = ind_start_index; i < N; i++)
    {
        const ChemicalSynapsisParams &params = individuals[i].params;

        synapsis.set(ChemicalSynapsisType::gfast, params.gfast);
        synapsis.set(ChemicalSynapsisType::gslow, params.gslow);
        synapsis.set(ChemicalSynapsisType::Esyn, params.Esyn);
        synapsis.set(ChemicalSynapsisType::sfast, params.sfast);
        synapsis.set(ChemicalSynapsisType::Vfast, params.Vfast);
        synapsis.set(ChemicalSynapsisType::Vslow, params.Vslow);
        synapsis.set(ChemicalSynapsisType::k1, params.k1);
        synapsis.set(ChemicalSynapsisType::k2, params.k2);
        synapsis.set(ChemicalSynapsisType::sslow, params.sslow);

        synapsis.set(ChemicalSynapsisType::mslow, FitnessConstants::M_SLOW_INITIAL_VALUE);

        reset_state_neur(model_neur);

        size_t interp_signal_counter = 0;
        size_t sig_i = 0;
        for (; sig_i < sig_start_minus_smoothed_avg_pts; sig_i++)
        {
            synapsis.step(dt, signal_data[sig_i], get_v_neur(model_neur));
            model_neur.add_synaptic_input(synapsis.get(i_enum));
            model_neur.step(dt);
            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_signal_data[interp_signal_counter], get_v_neur(model_neur));
                model_neur.add_synaptic_input(synapsis.get(i_enum));
                model_neur.step(dt);
                interp_signal_counter++;
            }
        }

        size_t model_sig_i = 0;
        for (; sig_i < signal_start_index; sig_i++, model_sig_i++)
        {
            const double v_neur = get_v_neur(model_neur);
            model_signal_buffer[model_sig_i] = v_neur;
            synapsis.step(dt, signal_data[sig_i], v_neur);
            model_neur.add_synaptic_input(synapsis.get(i_enum));
            model_neur.step(dt);
            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_signal_data[interp_signal_counter], get_v_neur(model_neur));
                model_neur.add_synaptic_input(synapsis.get(i_enum));
                model_neur.step(dt);
                interp_signal_counter++;
            }
        }

        size_t syn_sig_i = 0;
        for (; sig_i < total_size - 1; sig_i++, model_sig_i++, syn_sig_i++)
        {
            const double v_neur = get_v_neur(model_neur);
            model_signal_buffer[model_sig_i] = v_neur;
            synapsis.step(dt, signal_data[sig_i], v_neur);
            model_neur.add_synaptic_input(synapsis.get(i_enum));
            model_neur.step(dt);
            synapsis_signal_buffer[syn_sig_i] = synapsis.get(i_enum);
            if (use_ifast)
                ifast_signal_buffer[syn_sig_i] = synapsis.get(ifast_enum);
            if (use_islow)
                islow_signal_buffer[syn_sig_i] = synapsis.get(islow_enum);
            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_signal_data[interp_signal_counter], get_v_neur(model_neur));
                model_neur.add_synaptic_input(synapsis.get(i_enum));
                model_neur.step(dt);
                interp_signal_counter++;
            }
        }

        const double v_neur = get_v_neur(model_neur);
        model_signal_buffer[model_sig_i] = v_neur;
        synapsis.step(dt, signal_data[sig_i], v_neur);
        model_neur.add_synaptic_input(synapsis.get(i_enum));
        model_neur.step(dt);
        synapsis_signal_buffer[syn_sig_i] = synapsis.get(i_enum);
        if (use_ifast)
            ifast_signal_buffer[syn_sig_i] = synapsis.get(ifast_enum);
        if (use_islow)
            islow_signal_buffer[syn_sig_i] = synapsis.get(islow_enum);

        individuals[i].fitness = fitness_from_signals(living_const_signal_fitness_vals, search_phase, avg_smooth_points, use_ifast, use_islow, buffers);
    }
}