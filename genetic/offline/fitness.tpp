
#include <algorithm>
#include <cmath>
#include <iostream>
#include <kfr/all.hpp>
using namespace kfr;

namespace FitnessPrivateConfig
{
    static constexpr double SYN_RANGE_WEIGHT = 0.5;
    static constexpr double VPRE_SYN_COMP_WEIGHT = 0.5;

    static constexpr double IFAST_WEIGHT = 0.5;
    static constexpr double ISLOW_WEIGHT = 0.5;

    static constexpr double FILTER_FC = 0.3; // En KHz, es fija, la de las neuronas reales. Para que fs/fc=1000, cosa que vimos que es lo que debe salir, debemos poner un fs de 300 KHz, lo que implica un csv_step de 3.33... µs Usar 200000 puntos de esta duracion
    static constexpr double PAD_LEN_FACTOR = 1.5;
    static constexpr int BUTTERWORTH_ORDER = 4;
}

namespace FitnessConstants
{
    static constexpr double M_SLOW_INITIAL_VALUE = 0.0;
    static constexpr double I_MAX = 4.0;
    static constexpr double I_MIN = -4.0;
}

static double clamp_current_before_neuron_input(double i_val)
{
    constexpr double I_MAX = FitnessConstants::I_MAX;
    constexpr double I_MIN = FitnessConstants::I_MIN;
    return i_val >= I_MAX ? I_MAX : (i_val <= I_MIN ? I_MIN : i_val);
}

static double rescale_to_target(double value,
                                double src_min, double src_range,
                                double dst_min, double dst_range)
{
    const double norm = (value - src_min) / src_range;
    return dst_min + (norm * dst_range);
}

ConstantSigFitnessVals calc_const_sig_fitness_vals(
    const univector_ref<const double> &vpre_sig,
    double csv_step,
    bool use_ifast,
    bool use_islow,
    bool search_phase)
{
    ConstantSigFitnessVals result;

    const size_t use_size = vpre_sig.size();
    const double fs = 1000.0 / csv_step; // fs en KHz, csv_step en µs
    const double fc = FitnessPrivateConfig::FILTER_FC;

    const size_t effective_pad = std::min(static_cast<size_t>(FitnessPrivateConfig::PAD_LEN_FACTOR * fs / fc), use_size - 1);

    univector<double> padded(use_size + (2 * effective_pad));

    univector_ref<double> padded_seg = padded.slice(effective_pad, use_size);
    process(padded_seg, vpre_sig);

    double *padded_ptr = padded.data();
    const double *vpre_sig_ptr = vpre_sig.data();
    const double left_edge_x2 = 2.0 * vpre_sig_ptr[0];
    const double right_edge_x2 = 2.0 * vpre_sig_ptr[use_size - 1];
    for (size_t i = 0; i < effective_pad; i++)
    {
        padded_ptr[effective_pad - 1 - i] = left_edge_x2 - vpre_sig_ptr[i + 1];
        padded_ptr[use_size + effective_pad + i] = right_edge_x2 - vpre_sig_ptr[use_size - 2 - i];
    }

    filtfilt(padded, to_sos<double>(iir_lowpass(
                         butterworth(FitnessPrivateConfig::BUTTERWORTH_ORDER), fc, fs)));

    constexpr double EXPECTED_I_MIN_PHASE = FitnessPublicConfig::EXPECTED_I_MIN_PHASE;
    constexpr double EXPECTED_I_MAX_PHASE = FitnessPublicConfig::EXPECTED_I_MAX_PHASE;
    constexpr double EXPECTED_I_MIN_ANTIPHASE = FitnessPublicConfig::EXPECTED_I_MIN_ANTIPHASE;
    constexpr double EXPECTED_I_MAX_ANTIPHASE = FitnessPublicConfig::EXPECTED_I_MAX_ANTIPHASE;

    double expected_i_min;
    double expected_i_max;
    if (search_phase)
    {
        expected_i_min = EXPECTED_I_MIN_PHASE;
        expected_i_max = EXPECTED_I_MAX_PHASE;
    }
    else
    {
        expected_i_min = EXPECTED_I_MIN_ANTIPHASE;
        expected_i_max = EXPECTED_I_MAX_ANTIPHASE;
    }

    const bool use_both = use_ifast && use_islow;

    double vpre_min = 0.0, vpre_range = 0.0, expected_i_range = 0.0;
    if (use_both)
    {
        vpre_min = minof(vpre_sig);
        const double vpre_max = maxof(vpre_sig);
        vpre_range = vpre_max - vpre_min;
        expected_i_range = expected_i_max - expected_i_min;
    }

    if (use_islow)
    {
        univector<double> &islow_sig_centered = result.islow_sig_centered_to_fit;
        islow_sig_centered = padded_seg - mean(padded_seg);
        result.islow_sig_factor_to_fit = std::sqrt(sum(sqr(islow_sig_centered)));

        double &islow_sig_min = result.islow_sig_min_to_fit;
        double &islow_sig_max = result.islow_sig_max_to_fit;

        if (use_both)
        {
            islow_sig_min = rescale_to_target(minof(padded_seg), vpre_min, vpre_range,
                                              expected_i_min, expected_i_range);
            islow_sig_max = rescale_to_target(maxof(padded_seg), vpre_min, vpre_range,
                                              expected_i_min, expected_i_range);
        }
        else
        {
            islow_sig_min = expected_i_min;
            islow_sig_max = expected_i_max;
        }
        result.islow_sig_range_to_fit = islow_sig_max - islow_sig_min;
    }

    if (use_ifast)
    {
        univector<double> &ifast_sig_aux = result.ifast_sig_centered_to_fit;
        ifast_sig_aux = vpre_sig - padded_seg;

        double &ifast_sig_min = result.ifast_sig_min_to_fit;
        double &ifast_sig_max = result.ifast_sig_max_to_fit;

        if (use_both)
        {
            ifast_sig_min = rescale_to_target(minof(ifast_sig_aux), vpre_min, vpre_range,
                                              expected_i_min, expected_i_range);
            ifast_sig_max = rescale_to_target(maxof(ifast_sig_aux), vpre_min, vpre_range,
                                              expected_i_min, expected_i_range);
        }
        else
        {
            ifast_sig_min = expected_i_min;
            ifast_sig_max = expected_i_max;
        }
        result.ifast_sig_range_to_fit = ifast_sig_max - ifast_sig_min;

        ifast_sig_aux -= mean(ifast_sig_aux);
        result.ifast_sig_factor_to_fit = std::sqrt(sum(sqr(ifast_sig_aux)));
    }

    return result;
}

static double pearson_score(univector<double> &sig,
                            const univector<double> &ref_sig_centered,
                            double ref_sig_factor,
                            bool search_phase)
{
    sig -= mean(sig);
    const double sig_factor = std::sqrt(sum(sqr(sig)));
    const double r = sum(sig * ref_sig_centered) / (sig_factor * ref_sig_factor);
    const double normalized = (r + 1.0) / 2.0;
    return search_phase ? 1.0 - normalized : normalized;
}

static double calc_range_score_component(double observed_min, double observed_max,
                                         double expected_min, double expected_max,
                                         double max_diff)
{
    return 1.0 - (((std::abs(observed_min - expected_min) +
                    std::abs(observed_max - expected_max)) *
                   0.5) /
                  max_diff);
}

static inline double fitness_from_sigs(const ConstantSigFitnessVals &const_vpre_sig_fitness_vals, bool search_phase, bool use_ifast, bool use_islow, SigBuffers &buffers)
{
    constexpr double IFAST_WEIGHT = FitnessPrivateConfig::IFAST_WEIGHT;
    constexpr double ISLOW_WEIGHT = FitnessPrivateConfig::ISLOW_WEIGHT;

    double vpre_syn_comp_accum = 0.0;
    double i_range_accum = 0.0;
    double total_weight = 0.0;

    if (use_ifast)
    {
        univector<double> &ifast_sig = buffers.ifast_sig;
        const double comp_ifast = calc_range_score_component(
            minof(ifast_sig),
            maxof(ifast_sig),
            const_vpre_sig_fitness_vals.ifast_sig_min_to_fit,
            const_vpre_sig_fitness_vals.ifast_sig_max_to_fit,
            const_vpre_sig_fitness_vals.ifast_sig_range_to_fit);
        i_range_accum += IFAST_WEIGHT * comp_ifast;

        const double vpost_ifast_comp_score =
            pearson_score(ifast_sig,
                          const_vpre_sig_fitness_vals.ifast_sig_centered_to_fit,
                          const_vpre_sig_fitness_vals.ifast_sig_factor_to_fit,
                          search_phase);
        vpre_syn_comp_accum += IFAST_WEIGHT * vpost_ifast_comp_score;

        total_weight += IFAST_WEIGHT;
    }

    if (use_islow)
    {
        univector<double> &islow_sig = buffers.islow_sig;
        const double comp_islow = calc_range_score_component(
            minof(islow_sig),
            maxof(islow_sig),
            const_vpre_sig_fitness_vals.islow_sig_min_to_fit,
            const_vpre_sig_fitness_vals.islow_sig_max_to_fit,
            const_vpre_sig_fitness_vals.islow_sig_range_to_fit);
        i_range_accum += ISLOW_WEIGHT * comp_islow;

        const double vpost_islow_comp_score =
            pearson_score(islow_sig,
                          const_vpre_sig_fitness_vals.islow_sig_centered_to_fit,
                          const_vpre_sig_fitness_vals.islow_sig_factor_to_fit,
                          search_phase);
        vpre_syn_comp_accum += ISLOW_WEIGHT * vpost_islow_comp_score;

        total_weight += ISLOW_WEIGHT;
    }

    const double final_score = ((FitnessPrivateConfig::SYN_RANGE_WEIGHT * i_range_accum) +
                                (FitnessPrivateConfig::VPRE_SYN_COMP_WEIGHT * vpre_syn_comp_accum)) /
                               total_weight;

    if (!(std::isfinite(final_score)) || final_score < 0.0)
        return 0.0;

    return final_score;
}

template <typename Integrator, typename NeuronType, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    std::span<Individual> individuals,
                    const ScaledSigResult &scaled_result,
                    const ConstantSigFitnessVals &const_vpre_sig_fitness_vals,
                    bool search_phase,
                    SigBuffers &buffers,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur,
                    size_t vpre_sig_start_i,
                    bool use_ifast,
                    bool use_islow)
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    constexpr auto i_enum = ChemicalSynapsisType::i;
    constexpr auto ifast_enum = ChemicalSynapsisType::ifast;
    constexpr auto islow_enum = ChemicalSynapsisType::islow;

    double *ifast_sig_ptr = buffers.ifast_sig.data();
    double *islow_sig_ptr = buffers.islow_sig.data();

    const size_t total_size = scaled_result.sig.size();
    const size_t points_factor = scaled_result.points_factor;
    const double dt = scaled_result.dt;
    const double *vpre_sig_ptr = scaled_result.sig.data();
    const double *interpolated_points_ptr = scaled_result.interpolated_points.data();

    for (Individual &ind : individuals)
    {
        const ChemicalSynapsisVariationParams &params = ind.params;

        if (use_ifast)
        {
            synapsis.set(ChemicalSynapsisType::gfast, params.gfast);
            synapsis.set(ChemicalSynapsisType::sfast, params.sfast);
            synapsis.set(ChemicalSynapsisType::Vfast, params.Vfast);
        }

        if (use_islow)
        {
            synapsis.set(ChemicalSynapsisType::gslow, params.gslow);
            synapsis.set(ChemicalSynapsisType::Vslow, params.Vslow);
            synapsis.set(ChemicalSynapsisType::sslow, params.sslow);
            synapsis.set(ChemicalSynapsisType::k1, params.k1);
            synapsis.set(ChemicalSynapsisType::k2, params.k2);
        }

        synapsis.set(ChemicalSynapsisType::mslow, FitnessConstants::M_SLOW_INITIAL_VALUE);

        reset_state_neur(model_neur);

        size_t interp_pts_counter = 0;
        size_t vpre_sig_i = 0;
        for (; vpre_sig_i < vpre_sig_start_i; vpre_sig_i++)
        {
            synapsis.step(dt, vpre_sig_ptr[vpre_sig_i], get_v_neur(model_neur));
            const double i_val = clamp_current_before_neuron_input(synapsis.get(i_enum));
            model_neur.add_synaptic_input(-i_val);
            model_neur.step(dt);
            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_neur(model_neur));
                const double i_interp_val = clamp_current_before_neuron_input(synapsis.get(i_enum));
                model_neur.add_synaptic_input(-i_interp_val);
                model_neur.step(dt);
                interp_pts_counter++;
            }
        }

        size_t syn_sig_i = 0;
        for (; vpre_sig_i < total_size - 1; vpre_sig_i++, syn_sig_i++)
        {
            const double v_post = get_v_neur(model_neur);
            synapsis.step(dt, vpre_sig_ptr[vpre_sig_i], v_post);
            const double i_val = clamp_current_before_neuron_input(synapsis.get(i_enum));
            model_neur.add_synaptic_input(-i_val);
            model_neur.step(dt);
            if (use_ifast)
                ifast_sig_ptr[syn_sig_i] = synapsis.get(ifast_enum);
            if (use_islow)
                islow_sig_ptr[syn_sig_i] = synapsis.get(islow_enum);
            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_neur(model_neur));
                const double i_interp_val = clamp_current_before_neuron_input(synapsis.get(i_enum));
                model_neur.add_synaptic_input(-i_interp_val);
                model_neur.step(dt);
                interp_pts_counter++;
            }
        }

        const double v_post = get_v_neur(model_neur);
        synapsis.step(dt, vpre_sig_ptr[vpre_sig_i], v_post);
        const double i_val = clamp_current_before_neuron_input(synapsis.get(i_enum));
        model_neur.add_synaptic_input(-i_val);
        model_neur.step(dt);
        if (use_ifast)
            ifast_sig_ptr[syn_sig_i] = synapsis.get(ifast_enum);
        if (use_islow)
            islow_sig_ptr[syn_sig_i] = synapsis.get(islow_enum);

        ind.fitness = fitness_from_sigs(const_vpre_sig_fitness_vals, search_phase, use_ifast, use_islow, buffers);
    }
}