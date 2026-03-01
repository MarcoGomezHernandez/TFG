
#include <algorithm>
#include <cmath>
#include <iostream>
#include <kfr/all.hpp>
using namespace kfr;

namespace FitnessConfig
{
    static constexpr double V_COMP_WEIGHT = 0.5;
    static constexpr double VPRE_SYN_COMP_WEIGHT = 0.5;

    static constexpr double VPRE_I_COMP_WEIGHT = 0.4;
    static constexpr double VPRE_IFAST_COMP_WEIGHT = 0.3;
    static constexpr double VPRE_ISLOW_COMP_WEIGHT = 0.3;

    static constexpr double AVG_SMOOTH_POINTS_BURST_DIVISOR = 100;
    static constexpr double FILTER_FC = 24.76;

    static constexpr size_t FILTER_PAD_LEN = 1000;
    static constexpr int BUTTERWORTH_ORDER = 4;
}

namespace FitnessConstants
{
    static constexpr double M_SLOW_INITIAL_VALUE = 0.0;

    static constexpr size_t MIN_AVG_SMOOTH_POINTS = 1;
}

static void calc_syn_ref_vals(const univector_ref<const double> &vpre_sig,
                              ConstantSigFitnessVals &result,
                              double pts_burst_real,
                              bool use_ifast, bool use_islow,
                              double min, double max)
{
    const size_t use_size = vpre_sig.size();
    const double fs = pts_burst_real;
    const double fc = FitnessConfig::FILTER_FC;

    const size_t effective_pad = std::min(FitnessConfig::FILTER_PAD_LEN, use_size - 1);

    univector<double> padded(use_size + (2 * effective_pad));

    univector_ref<double> padded_seg = padded.slice(effective_pad, use_size);
    process(padded_seg, vpre_sig);

    double *padded_ptr = padded.data();
    const double *vpre_sig_ptr = vpre_sig.data();
    for (size_t i = 0; i < effective_pad; i++)
    {
        padded_ptr[effective_pad - 1 - i] = vpre_sig_ptr[i + 1];
        padded_ptr[use_size + effective_pad + i] = vpre_sig_ptr[use_size - 2 - i];
    }

    filtfilt(padded, to_sos<double>(iir_lowpass(
                         butterworth(FitnessConfig::BUTTERWORTH_ORDER), fc, fs)));

    if (use_ifast && use_islow)
    {
        univector<double> &i_sig_centered = result.i_sig_centered_to_fit;
        i_sig_centered = vpre_sig - mean(vpre_sig);
        result.i_sig_factor_to_fit = std::sqrt(sum(sqr(i_sig_centered)));
    }

    if (use_islow)
    {
        univector<double> &islow_sig_centered = result.islow_sig_centered_to_fit;
        islow_sig_centered = padded_seg - mean(padded_seg);
        result.islow_sig_factor_to_fit = std::sqrt(sum(sqr(islow_sig_centered)));
    }

    if (use_ifast)
    {
        univector<double> ifast_sig = vpre_sig - padded_seg;
        univector<double> &ifast_sig_centered = result.ifast_sig_centered_to_fit;
        ifast_sig_centered = ifast_sig - mean(ifast_sig);
        result.ifast_sig_factor_to_fit = std::sqrt(sum(sqr(ifast_sig_centered)));
    }
}

ConstantSigFitnessVals calc_const_sig_fitness_vals(
    const univector<double> &vpre_sig,
    double min,
    double max,
    bool search_phase,
    size_t avg_smooth_points,
    double pts_burst_real,
    bool use_ifast,
    bool use_islow)
{
    ConstantSigFitnessVals result;

    const size_t vpre_sig_size = vpre_sig.size();
    const size_t use_size = vpre_sig_size - avg_smooth_points;

    univector<double> &smoothed_vpre_sig_to_fit = result.smoothed_vpre_sig_to_fit;
    smoothed_vpre_sig_to_fit.reserve(use_size);

    double running_sum = sum(vpre_sig.slice(0, avg_smooth_points));

    const double *vpre_sig_ptr = vpre_sig.data();
    for (size_t i = avg_smooth_points; i < vpre_sig_size; i++)
    {
        running_sum += vpre_sig_ptr[i];
        running_sum -= vpre_sig_ptr[i - avg_smooth_points];
        smoothed_vpre_sig_to_fit.push_back(running_sum / avg_smooth_points);
    }

    const double smoothed_min = minof(smoothed_vpre_sig_to_fit);
    const double smoothed_max = maxof(smoothed_vpre_sig_to_fit);

    if (!search_phase)
    {
        smoothed_vpre_sig_to_fit = (smoothed_min + smoothed_max) - smoothed_vpre_sig_to_fit;
    }

    result.max_v_comp_distance = smoothed_max - smoothed_min;

    calc_syn_ref_vals(vpre_sig.slice(avg_smooth_points, use_size),
                      result,
                      pts_burst_real,
                      use_ifast,
                      use_islow,
                      min,
                      max);

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
    return search_phase ? normalized : 1.0 - normalized;
}

static double fitness_from_sigs(const ConstantSigFitnessVals &const_vpre_sig_fitness_vals, bool search_phase, size_t avg_smooth_points, bool use_ifast, bool use_islow, SigBuffers &buffers)
{
    const univector<double> &vpost_sig = buffers.vpost_sig;
    univector<double> &smoothed_vpost_sig = buffers.smoothed_vpost_sig;
    const size_t vpost_sig_size = vpost_sig.size();

    double running_sum = sum(vpost_sig.slice(0, avg_smooth_points));

    const double *vpost_sig_ptr = vpost_sig.data();
    double *smoothed_vpost_sig_ptr = smoothed_vpost_sig.data();
    for (size_t i = avg_smooth_points; i < vpost_sig_size; i++)
    {
        running_sum += vpost_sig_ptr[i];
        running_sum -= vpost_sig_ptr[i - avg_smooth_points];
        smoothed_vpost_sig_ptr[i - avg_smooth_points] = running_sum / avg_smooth_points;
    }
    const double v_comp_RMSE = std::sqrt(sum(sqr(smoothed_vpost_sig - const_vpre_sig_fitness_vals.smoothed_vpre_sig_to_fit)) / smoothed_vpost_sig.size());
    const double v_comp_score = 1.0 - (v_comp_RMSE / const_vpre_sig_fitness_vals.max_v_comp_distance);

    double vpre_syn_comp_score = 0.0;
    double vpre_syn_comp_weight = 0.0;

    if (use_ifast && use_islow)
    {
        const double vpre_i_comp_score =
            pearson_score(buffers.i_sig,
                          const_vpre_sig_fitness_vals.i_sig_centered_to_fit,
                          const_vpre_sig_fitness_vals.i_sig_factor_to_fit,
                          search_phase);
        vpre_syn_comp_score += FitnessConfig::VPRE_I_COMP_WEIGHT * vpre_i_comp_score;
        vpre_syn_comp_weight += FitnessConfig::VPRE_I_COMP_WEIGHT;
    }

    if (use_ifast)
    {
        const double vpost_ifast_comp_score =
            pearson_score(buffers.ifast_sig,
                          const_vpre_sig_fitness_vals.ifast_sig_centered_to_fit,
                          const_vpre_sig_fitness_vals.ifast_sig_factor_to_fit,
                          search_phase);
        vpre_syn_comp_score += FitnessConfig::VPRE_IFAST_COMP_WEIGHT * vpost_ifast_comp_score;
        vpre_syn_comp_weight += FitnessConfig::VPRE_IFAST_COMP_WEIGHT;
    }

    if (use_islow)
    {
        const double vpost_islow_comp_score =
            pearson_score(buffers.islow_sig,
                          const_vpre_sig_fitness_vals.islow_sig_centered_to_fit,
                          const_vpre_sig_fitness_vals.islow_sig_factor_to_fit,
                          search_phase);
        vpre_syn_comp_score += FitnessConfig::VPRE_ISLOW_COMP_WEIGHT * vpost_islow_comp_score;
        vpre_syn_comp_weight += FitnessConfig::VPRE_ISLOW_COMP_WEIGHT;
    }

    const double final_score = (FitnessConfig::V_COMP_WEIGHT * v_comp_score) +
                               (FitnessConfig::VPRE_SYN_COMP_WEIGHT * (vpre_syn_comp_score / vpre_syn_comp_weight));

    if (!(std::isfinite(final_score)) || final_score < 0.0)
        return 0.0;

    return final_score;
}

template <typename Integrator, typename NeuronType, size_t N, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    std::array<Individual, N> &individuals,
                    const ScaledSigResult &scaled_result,
                    const ConstantSigFitnessVals &const_vpre_sig_fitness_vals,
                    bool search_phase,
                    SigBuffers &buffers,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur,
                    size_t ind_start_i,
                    size_t vpre_sig_start_i,
                    size_t avg_smooth_points,
                    bool use_ifast,
                    bool use_islow)
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    constexpr auto i_enum = ChemicalSynapsisType::i;
    constexpr auto ifast_enum = ChemicalSynapsisType::ifast;
    constexpr auto islow_enum = ChemicalSynapsisType::islow;

    double *vpost_sig_ptr = buffers.vpost_sig.data();
    double *i_sig_ptr = buffers.i_sig.data();
    double *ifast_sig_ptr = buffers.ifast_sig.data();
    double *islow_sig_ptr = buffers.islow_sig.data();

    const size_t total_size = scaled_result.sig.size();
    const size_t points_factor = scaled_result.points_factor;
    const double dt = scaled_result.dt;
    const double *vpre_sig_ptr = scaled_result.sig.data();
    const double *interpolated_points_ptr = scaled_result.interpolated_points.data();

    for (size_t i = ind_start_i; i < N; i++)
    {
        Individual &ind = individuals[i];
        const ChemicalSynapsisParams &params = ind.params;

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

        size_t interp_pts_counter = 0;
        size_t vpre_sig_i = 0;
        for (; vpre_sig_i < vpre_sig_start_i - avg_smooth_points; vpre_sig_i++)
        {
            synapsis.step(dt, vpre_sig_ptr[vpre_sig_i], get_v_neur(model_neur));
            model_neur.add_synaptic_input(synapsis.get(i_enum));
            model_neur.step(dt);
            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_neur(model_neur));
                model_neur.add_synaptic_input(synapsis.get(i_enum));
                model_neur.step(dt);
                interp_pts_counter++;
            }
        }

        size_t vpost_sig_i = 0;
        for (; vpre_sig_i < vpre_sig_start_i; vpre_sig_i++, vpost_sig_i++)
        {
            const double v_post = get_v_neur(model_neur);
            vpost_sig_ptr[vpost_sig_i] = v_post;
            synapsis.step(dt, vpre_sig_ptr[vpre_sig_i], v_post);
            model_neur.add_synaptic_input(synapsis.get(i_enum));
            model_neur.step(dt);
            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_neur(model_neur));
                model_neur.add_synaptic_input(synapsis.get(i_enum));
                model_neur.step(dt);
                interp_pts_counter++;
            }
        }

        size_t syn_sig_i = 0;
        for (; vpre_sig_i < total_size - 1; vpre_sig_i++, vpost_sig_i++, syn_sig_i++)
        {
            const double v_post = get_v_neur(model_neur);
            vpost_sig_ptr[vpost_sig_i] = v_post;
            synapsis.step(dt, vpre_sig_ptr[vpre_sig_i], v_post);
            const double i_val = synapsis.get(i_enum);
            model_neur.add_synaptic_input(i_val);
            model_neur.step(dt);
            if (use_ifast && use_islow)
                i_sig_ptr[syn_sig_i] = i_val;
            if (use_ifast)
                ifast_sig_ptr[syn_sig_i] = synapsis.get(ifast_enum);
            if (use_islow)
                islow_sig_ptr[syn_sig_i] = synapsis.get(islow_enum);
            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_neur(model_neur));
                model_neur.add_synaptic_input(synapsis.get(i_enum));
                model_neur.step(dt);
                interp_pts_counter++;
            }
        }

        const double v_post = get_v_neur(model_neur);
        vpost_sig_ptr[vpost_sig_i] = v_post;
        synapsis.step(dt, vpre_sig_ptr[vpre_sig_i], v_post);
        const double i_val = synapsis.get(i_enum);
        model_neur.add_synaptic_input(i_val);
        model_neur.step(dt);
        if (use_ifast && use_islow)
            i_sig_ptr[syn_sig_i] = i_val;
        if (use_ifast)
            ifast_sig_ptr[syn_sig_i] = synapsis.get(ifast_enum);
        if (use_islow)
            islow_sig_ptr[syn_sig_i] = synapsis.get(islow_enum);

        ind.fitness = fitness_from_sigs(const_vpre_sig_fitness_vals, search_phase, avg_smooth_points, use_ifast, use_islow, buffers);
    }
}