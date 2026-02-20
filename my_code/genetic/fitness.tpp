
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
 * start_index: leading points used only for rolling-average warm-up (not stored in output).
 *   Must satisfy start_index >= avg_smooth_points so the window is already full at start_index,
 *   allowing the recorded portion to start directly in steady-state (no warm-up sub-loop needed).
 * Pass 1: causal running-average smooth + compute min/max of smoothed signal.
 * Pass 2 (antiphase only): flip smoothed signal around [min, max].
 * Pass 3: normalize.
 */
ConstantSignalFitnessVals calc_const_signal_fitness_vals(const std::vector<double> &signal, double min_val, double max_val, bool search_phase, size_t avg_smooth_points, size_t start_index)
{
    ConstantSignalFitnessVals result;

    const size_t total_size = signal.size();
    const size_t use_size = total_size - start_index;

    std::vector<double> &smoothed_signal_to_fit = result.smoothed_signal_to_fit;
    smoothed_signal_to_fit.resize(use_size);
    std::vector<double> &normalized_signal_to_fit = result.normalized_signal_to_fit;
    normalized_signal_to_fit.resize(use_size);

    // Pass 1: run rolling-average over stabilization prefix [0..start_index-1] without storing,
    //         then compute smoothed+normalized over use portion [start_index..total_size-1].
    // Since start_index >= avg_smooth_points, the window is full by start_index, so no warm-up
    // sub-loop is needed inside the recorded portion.
    const double range = max_val - min_val;
    double running_sum = 0.0;
    double smoothed_min_val = GeneralConstants::DOUBLE_MAX;
    double smoothed_max_val = GeneralConstants::DOUBLE_MIN;

    // Pre-fill window: sum the avg_smooth_points values immediately before start_index.
    // (Guaranteed valid since start_index >= avg_smooth_points.)
    for (size_t i = start_index - avg_smooth_points; i < start_index; i++)
        running_sum += signal[i];

    // Main loop: record use portion [start_index..total_size-1] in steady state (no warm-up)
    for (size_t i = start_index; i < total_size; i++)
    {
        const double sig_val = signal[i];

        running_sum += sig_val;
        running_sum -= signal[i - avg_smooth_points];
        const double smoothed_val = running_sum / avg_smooth_points;
        const size_t i_0_start = i - start_index;
        smoothed_signal_to_fit[i_0_start] = smoothed_val;

        normalized_signal_to_fit[i_0_start] = (sig_val - min_val) / range;

        if (smoothed_val < smoothed_min_val)
            smoothed_min_val = smoothed_val;
        if (smoothed_val > smoothed_max_val)
            smoothed_max_val = smoothed_val;
    }

    // Precompute maximum possible v_comp distance for normalization:
    //   max_v_comp_distance = use_size * (max_smoothed_val - min_smoothed_val)
    result.max_v_comp_distance = use_size * (smoothed_max_val - smoothed_min_val);

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
 * model_signal layout: [avg_smooth_points_model seed values | use_size model voltages]
 *   The seed is the last avg_smooth_points_model outer-loop model voltages from the stabilization
 *   phase (chronologically ordered by calc_fitnesses). This pre-fills the rolling-average window
 *   so that the comparison loop runs entirely in steady state (no warm-up sub-loop needed).
 * Returns the computed fitness score
 */
double fitness_from_signals(const ConstantSignalFitnessVals &living_const_signal_fitness_vals, const std::vector<double> &model_signal, bool search_phase, const std::vector<double> &synapsis_signal, size_t avg_smooth_points_model)
{
    // Initialize running_sum from the seed prefix (chronologically ordered stab tail)
    double running_sum = 0.0;
    for (size_t i = 0; i < avg_smooth_points_model; i++)
        running_sum += model_signal[i];

    const size_t use_size = synapsis_signal.size();

    // Compute v_comp_score on-the-fly: steady loop, no warm-up needed
    double v_comp_dist = 0.0;
    const std::vector<double> &living_smoothed_signal_to_fit = living_const_signal_fitness_vals.smoothed_signal_to_fit;
    for (size_t i = 0; i < use_size; i++)
    {
        running_sum += model_signal[i + avg_smooth_points_model];
        running_sum -= model_signal[i];
        const double smoothed_val = running_sum / avg_smooth_points_model;
        v_comp_dist += std::abs(living_smoothed_signal_to_fit[i] - smoothed_val);
    }
    const double v_comp_score = 1.0 - (v_comp_dist / living_const_signal_fitness_vals.max_v_comp_distance); // Normalize: max possible distance is (max-min)*use_size

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
    for (size_t i = 0; i < use_size; i++)
    {
        const double norm_syn_val = (synapsis_signal[i] - syn_min) / syn_range;
        vpre_i_comp_dist += std::abs(living_norm_signal_to_fit[i] - norm_syn_val);
    }
    // Normalize: max possible accumulated distance is 1.0 * use_size (both signals in [0,1])
    const double vpre_i_comp_score = 1.0 - (vpre_i_comp_dist / use_size);

    // Compute final weighted score
    const double final_score = (FitnessConfig::V_COMP_WEIGHT * v_comp_score) + (FitnessConfig::VPRE_I_COMP_WEIGHT * vpre_i_comp_score);

    // If the computed fitness is NaN, return 0.0 instead
    if (std::isnan(final_score) || final_score < 0.0)
        return 0.0;

    return final_score;
}

/*
 * Template function to calculate fitnesses for multiple parameter sets
 * Simulates the neural model with given parameters and computes fitness against precomputed vals
 * Parameters: synapsis, neurons, individuals, scaled_result, living_const_signal_fitness_vals, search_phase, buffers, reset_state_neur, get_v_neur, ind_start_index, signal_start_index, avg_smooth_points_model
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
                    size_t ind_start_index,
                    size_t signal_start_index,
                    size_t avg_smooth_points_model)
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    const size_t signal_size = scaled_result.signal.size();
    const size_t points_factor = scaled_result.points_factor;
    const double dt = scaled_result.dt;
    const double *signal_data = scaled_result.signal.data();
    const double *interpolated_signal_data = scaled_result.interpolated_points.data();

    const size_t sig_start_minus_smoothed_avg_pts = signal_start_index - avg_smooth_points_model;

    for (size_t i = ind_start_index; i < N; i++)
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

        for (; j < sig_start_minus_smoothed_avg_pts; j++)
        {
            synapsis.step(dt, signal_data[j], get_v_neur(model_neur));
            model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
            model_neur.step(dt);
            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_signal_data[interp_signal_counter], get_v_neur(model_neur));
                model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
                model_neur.step(dt);
                interp_signal_counter++;
            }
        }

        for (; j < signal_start_index; j++)
        {
            const double v_neur = get_v_neur(model_neur);
            model_signal_buffer[j] = v_neur;
            synapsis.step(dt, signal_data[j], v_neur);
            model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
            model_neur.step(dt);
            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_signal_data[interp_signal_counter], get_v_neur(model_neur));
                model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
                model_neur.step(dt);
                interp_signal_counter++;
            }
        }

        for (; j < signal_size - 1; j++)
        {
            const double v_neur = get_v_neur(model_neur);
            model_signal_buffer[j] = v_neur;
            synapsis.step(dt, signal_data[j], v_neur);
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

        const double v_neur = get_v_neur(model_neur);
        model_signal_buffer[j] = v_neur;
        synapsis.step(dt, signal_data[j], v_neur);
        model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
        model_neur.step(dt);
        synapsis_signal_buffer[j] = synapsis.get(ChemicalSynapsisType::i);

        // Compute fitness and store directly in individual
        individuals[i].fitness = fitness_from_signals(living_const_signal_fitness_vals, model_signal_buffer, search_phase, synapsis_signal_buffer, avg_smooth_points_model);
    }
}