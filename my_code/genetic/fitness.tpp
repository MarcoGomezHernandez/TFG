
#include <algorithm>
#include <cmath>
#include <iostream>

/*
 * Configuration constants for fitness calculation
 */
namespace FitnessConfig
{
    // Weights for score components (must sum to 1.0)
    static constexpr double V_COMP_WEIGHT = 0.5;
    static constexpr double VPRE_I_COMP_WEIGHT = 0.5;

    // Fraction of points-per-burst used as the smoothing window
    static constexpr double AVG_SMOOTH_POINTS_BURST_FRACTION = 0.01;
}

namespace FitnessConstants
{
    static constexpr double M_SLOW_INITIAL_VALUE = 0.0;

    // Minimum smoothing window size (no smoothing below this)
    static constexpr size_t MIN_AVG_SMOOTH_POINTS = 1;
}

/*
 * Function to preprocess a signal and compute statistics
 * Returns ConstantSignalFitnessVals containing normalized_signal_to_fit and smoothed_signal_to_fit
 * search_phase: true for v_comp (signal as-is), false for antiphase (signal flipped around min/max)
 * avg_smooth_points: number of points to average for smoothing (1 = no smoothing)
 * Pass 1: causal running-average smooth + compute min/max of smoothed signal.
 * Pass 2 (antiphase only): flip smoothed signal around [min, max].
 * Pass 3: normalize.
 */
ConstantSignalFitnessVals calc_const_signal_fitness_vals(const std::vector<double> &signal, double min_val, double max_val, bool search_phase, size_t avg_smooth_points)
{
    ConstantSignalFitnessVals result;

    const size_t signal_size = signal.size();

    std::vector<double> &smoothed_signal_to_fit = result.smoothed_signal_to_fit;
    smoothed_signal_to_fit.resize(signal_size);
    std::vector<double> &normalized_signal_to_fit = result.normalized_signal_to_fit;
    normalized_signal_to_fit.resize(signal_size);

    // Pass 1: causal running-average smooth + normalize using known [min_val, max_val]
    const double range = max_val - min_val;
    double running_sum = 0.0;
    double smoothed_min_val = GeneralConstants::DOUBLE_MAX;
    double smoothed_max_val = GeneralConstants::DOUBLE_MIN;
    // First sub-loop: warm-up v_comp (window grows, no subtraction needed)
    const size_t warm_up = std::min(avg_smooth_points, signal_size);
    for (size_t i = 0; i < warm_up; i++)
    {
        running_sum += signal[i];
        const double smoothed_val = running_sum / (i + 1);
        smoothed_signal_to_fit[i] = smoothed_val;
        normalized_signal_to_fit[i] = (signal[i] - min_val) / range;
        if (smoothed_val < smoothed_min_val)
            smoothed_min_val = smoothed_val;
        if (smoothed_val > smoothed_max_val)
            smoothed_max_val = smoothed_val;
    }
    // Second sub-loop: steady v_comp (window = avg_smooth_points, subtract oldest)
    for (size_t i = warm_up; i < signal_size; i++)
    {
        running_sum += signal[i];
        running_sum -= smoothed_signal_to_fit[i - avg_smooth_points];
        const double smoothed_val = running_sum / avg_smooth_points;
        smoothed_signal_to_fit[i] = smoothed_val;
        normalized_signal_to_fit[i] = (signal[i] - min_val) / range;
        if (smoothed_val < smoothed_min_val)
            smoothed_min_val = smoothed_val;
        if (smoothed_val > smoothed_max_val)
            smoothed_max_val = smoothed_val;
    }

    // Precompute maximum possible v_comp distance for normalization:
    //   max_v_comp_distance = signal_size * (max_smoothed_val - min_smoothed_val)
    result.max_v_comp_distance = signal_size * (smoothed_max_val - smoothed_min_val);

    // Pass 2 (antiphase only): flip both signals
    if (!search_phase)
    {
        // Flip smoothed signal around its own [smoothed_min_val, smoothed_max_val]
        const double smoothed_flip_offset = smoothed_min_val + smoothed_max_val;
        for (double &val : smoothed_signal_to_fit)
            val = smoothed_flip_offset - val;

        // Flip normalized signal around [min_val, max_val] → equivalent to 1 - val
        for (double &val : normalized_signal_to_fit)
            val = 1.0 - val;
    }

    return result;
}

/*
 * Calculate fitness as distance-based scores between model/synapsis signals and the reference signal
 * Takes precomputed vals for the reference signal, raw model_signal, synapsis_signal
 * avg_smooth_points_model: smoothing window for the model signal (1 = no smoothing)
 * Returns the computed fitness score
 */
double fitness_from_signals(const ConstantSignalFitnessVals &living_const_signal_fitness_vals, const std::vector<double> &model_signal, bool search_phase, const std::vector<double> &synapsis_signal, size_t avg_smooth_points_model)
{
    const size_t signal_size = model_signal.size();

    // Compute v_comp_score on-the-fly using a causal running average of the model signal
    // First sub-loop: warm-up v_comp (window grows, no subtraction needed)
    double v_comp_dist = 0.0;
    double running_sum = 0.0;
    const std::vector<double> &living_smoothed_signal_to_fit = living_const_signal_fitness_vals.smoothed_signal_to_fit;
    const size_t warm_up = std::min(avg_smooth_points_model, signal_size);
    for (size_t i = 0; i < warm_up; i++)
    {
        running_sum += model_signal[i];
        const double smoothed_val = running_sum / (i + 1);
        v_comp_dist += std::abs(living_smoothed_signal_to_fit[i] - smoothed_val);
    }
    // Second sub-loop: steady v_comp (window = avg_smooth_points_model, subtract oldest)
    for (size_t i = warm_up; i < signal_size; i++)
    {
        running_sum += model_signal[i];
        running_sum -= model_signal[i - avg_smooth_points_model];
        const double smoothed_val = running_sum / avg_smooth_points_model;
        v_comp_dist += std::abs(living_smoothed_signal_to_fit[i] - smoothed_val);
    }
    const double v_comp_score = 1.0 - (v_comp_dist / living_const_signal_fitness_vals.max_v_comp_distance); // Normalize: max possible distance is (max-min)*signal_size

    // Compute vpre_i_comp_score: normalized distance between synapsis signal and reference signal_to_fit
    double syn_min = GeneralConstants::DOUBLE_MAX;
    double syn_max = GeneralConstants::DOUBLE_MIN;
    for (double syn_val : synapsis_signal)
    {
        if (syn_val < syn_min)
            syn_min = syn_val;
        if (syn_val > syn_max)
            syn_max = syn_val;
    }

    const std::vector<double> &living_norm_signal_to_fit = living_const_signal_fitness_vals.normalized_signal_to_fit;
    const double syn_range = syn_max - syn_min;
    double vpre_i_comp_dist = 0.0;
    for (size_t i = 0; i < signal_size; i++)
    {
        const double norm_syn_val = (synapsis_signal[i] - syn_min) / syn_range;
        vpre_i_comp_dist += std::abs(living_norm_signal_to_fit[i] - norm_syn_val);
    }
    // Normalize: max possible accumulated distance is 1.0 * signal_size (both signals in [0,1])
    const double vpre_i_comp_score = 1.0 - (vpre_i_comp_dist / signal_size);

    // Compute final weighted score
    const double final_score = (FitnessConfig::V_COMP_WEIGHT * v_comp_score) + (FitnessConfig::VPRE_I_COMP_WEIGHT * vpre_i_comp_score);

    // If the computed fitness is NaN, return 0.0 instead
    if (final_score < 0.0)
        return 0.0;

    return final_score;
}

/*
 * Template function to calculate fitnesses for multiple parameter sets
 * Simulates the neural model with given parameters and computes fitness against precomputed vals
 * Parameters: synapsis, neurons, individuals, scaled_result, living_const_signal_fitness_vals, search_phase, buffers, reset_state_neur, get_v_neur, start_index, avg_smooth_points_model
 */
template <typename Integrator, typename NeuronType, size_t N, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    std::array<Individual, N> &individuals,
                    const ScaledSignalResult &scaled_result,
                    const ConstantSignalFitnessVals &living_const_signal_fitness_vals,
                    bool search_phase,
                    std::vector<double> &model_signal_buffer,
                    std::vector<double> &synapsis_signal_buffer,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur,
                    size_t start_index,
                    size_t avg_smooth_points_model)
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    const size_t signal_size = scaled_result.signal.size();
    const size_t points_factor = scaled_result.points_factor;
    const double dt = scaled_result.dt;
    const double *signal_data = scaled_result.signal.data();
    const double *interpolated_signal_data = scaled_result.interpolated_points.data();

    for (size_t i = start_index; i < N; i++)
    {
        const ChemicalSynapsisParams &params = individuals[i].params;

        // Set synapsis parameters
        synapsis.set(ChemicalSynapsisType::gfast, params.gfast);
        synapsis.set(ChemicalSynapsisType::gslow, params.gslow);
        synapsis.set(ChemicalSynapsisType::Esyn, params.Esyn);
        synapsis.set(ChemicalSynapsisType::sfast, params.sfast);
        synapsis.set(ChemicalSynapsisType::Vfast, params.Vfast);
        synapsis.set(ChemicalSynapsisType::Vslow, params.Vslow);
        synapsis.set(ChemicalSynapsisType::k1, params.k1);
        synapsis.set(ChemicalSynapsisType::k2, params.k2);
        synapsis.set(ChemicalSynapsisType::sslow, params.sslow);

        // Reset synapsis
        synapsis.set(ChemicalSynapsisType::mslow, FitnessConstants::M_SLOW_INITIAL_VALUE);

        // Reset neuron state using provided function
        reset_state_neur(model_neur);

        // Simulate and collect neur signal and synapsis current
        size_t interp_signal_counter = 0;
        size_t j = 0;
        for (; j < signal_size - 1; j++)
        {
            model_signal_buffer[j] = get_v_neur(model_neur);
            synapsis.step(dt, signal_data[j], get_v_neur(model_neur));
            model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
            model_neur.step(dt);
            synapsis_signal_buffer[j] = synapsis.get(ChemicalSynapsisType::i);

            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_signal_data[interp_signal_counter], get_v_neur(model_neur));
                model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
                model_neur.step(dt);
                interp_signal_counter++;
            }
        }

        model_signal_buffer[j] = get_v_neur(model_neur);
        synapsis.step(dt, signal_data[j], get_v_neur(model_neur));
        model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
        model_neur.step(dt);
        synapsis_signal_buffer[j] = synapsis.get(ChemicalSynapsisType::i);

        // Compute fitness and store directly in individual
        individuals[i].fitness = fitness_from_signals(living_const_signal_fitness_vals, model_signal_buffer, search_phase, synapsis_signal_buffer, avg_smooth_points_model);
    }
}