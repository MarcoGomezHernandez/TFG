
#include <algorithm>
#include <cmath>
#include <iostream>
#include <kfr/all.hpp>
using namespace kfr;

namespace FitnessConfig
{
    static constexpr double SYN_RANGE_WEIGHT = 0.5;
    static constexpr double VPRE_SYN_COMP_WEIGHT = 0.5;

    static constexpr double VPRE_I_COMP_WEIGHT = 0.4;
    static constexpr double VPRE_IFAST_COMP_WEIGHT = 0.3;
    static constexpr double VPRE_ISLOW_COMP_WEIGHT = 0.3;

    static constexpr double SYN_RANGE_EXPECTED_MIN_PHASE = 0.25;
    static constexpr double SYN_RANGE_EXPECTED_MAX_PHASE = 1.1;

    static constexpr double SYN_RANGE_EXPECTED_MIN_ANTIPHASE = 0.25;
    static constexpr double SYN_RANGE_EXPECTED_MAX_ANTIPHASE = 1.1;

    static constexpr double FILTER_FC = 24.76;
    static constexpr size_t FILTER_PAD_LEN = 1000;
    static constexpr int BUTTERWORTH_ORDER = 4;
}

namespace FitnessConstants
{
    static constexpr double M_SLOW_INITIAL_VALUE = 0.0;

    static constexpr double SYN_RANGE_MAX_DIFF_PHASE = FitnessConfig::SYN_RANGE_EXPECTED_MAX_PHASE - FitnessConfig::SYN_RANGE_EXPECTED_MIN_PHASE;
    static constexpr double SYN_RANGE_MAX_DIFF_ANTIPHASE = FitnessConfig::SYN_RANGE_EXPECTED_MAX_ANTIPHASE - FitnessConfig::SYN_RANGE_EXPECTED_MIN_ANTIPHASE;
}

ConstantSigFitnessVals calc_const_sig_fitness_vals(
    const univector_ref<const double> &vpre_sig,
    double pts_burst_real,
    bool use_ifast,
    bool use_islow)
{
    ConstantSigFitnessVals result;

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

    return result;
}

static inline double pearson_score(univector<double> &sig,
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

static inline double range_score(const SigBuffers &buffers, bool use_ifast, bool use_islow, bool search_phase)
{
    double observed_min;
    double observed_max;

    if (use_ifast && use_islow)
    {
        observed_min = minof(buffers.i_sig);
        observed_max = maxof(buffers.i_sig);
    }
    else if (use_ifast)
    {
        observed_min = minof(buffers.ifast_sig);
        observed_max = maxof(buffers.ifast_sig);
    }
    else
    {
        observed_min = minof(buffers.islow_sig);
        observed_max = maxof(buffers.islow_sig);
    }

    double expected_min, expected_max, max_diff;
    if (search_phase)
    {
        expected_min = FitnessConfig::SYN_RANGE_EXPECTED_MIN_PHASE;
        expected_max = FitnessConfig::SYN_RANGE_EXPECTED_MAX_PHASE;
        max_diff = FitnessConstants::SYN_RANGE_MAX_DIFF_PHASE;
    }
    else
    {
        expected_min = FitnessConfig::SYN_RANGE_EXPECTED_MIN_ANTIPHASE;
        expected_max = FitnessConfig::SYN_RANGE_EXPECTED_MAX_ANTIPHASE;
        max_diff = FitnessConstants::SYN_RANGE_MAX_DIFF_ANTIPHASE;
    }

    return 1.0 - (((std::abs(observed_min - expected_min) +
                    std::abs(observed_max - expected_max)) *
                   0.5) /
                  max_diff);
}

static inline double fitness_from_sigs(const ConstantSigFitnessVals &const_vpre_sig_fitness_vals, bool search_phase, bool use_ifast, bool use_islow, SigBuffers &buffers)
{
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

    const double syn_range_score = range_score(buffers, use_ifast, use_islow, search_phase);

    const double final_score = (FitnessConfig::SYN_RANGE_WEIGHT * syn_range_score) +
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
                    bool use_ifast,
                    bool use_islow)
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    constexpr auto i_enum = ChemicalSynapsisType::i;
    constexpr auto ifast_enum = ChemicalSynapsisType::ifast;
    constexpr auto islow_enum = ChemicalSynapsisType::islow;

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
        for (; vpre_sig_i < total_size - 1; vpre_sig_i++, syn_sig_i++)
        {
            const double v_post = get_v_neur(model_neur);
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

        ind.fitness = fitness_from_sigs(const_vpre_sig_fitness_vals, search_phase, use_ifast, use_islow, buffers);
    }
}