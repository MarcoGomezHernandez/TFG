#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "utils.h"
#include "scaling.h"
#include <ChemicalSynapsis.h>

/*
 * Configuration constants for fitness calculation
 */
namespace FitnessConfig
{
    // Weights for score components
    static constexpr double BURSTS_DIFF_WEIGHT = 0.3;
    static constexpr double PHASE_WEIGHT = 0.4;
    static constexpr double MINMAX_WEIGHT = 0.3;
}

namespace FitnessConstants
{
    static constexpr double NORM_MAX_MINMAX_DIFF = (HindmarshRose::MAX - HindmarshRose::MIN) * 2.0;
    static constexpr double M_SLOW_INITIAL_VALUE = 0.0;
}

/*
 * Compute an inverse normalization that maps a value to a score between 0 and 1, with score=1 at val=min_val and score=0 at val=max_val, clamped to [0,1].
 */
double inverse_normalization(double val, double min_val, double max_val)
{
    return (max_val - val) / (max_val - min_val);
}

/*
 * Struct to hold precomputed signal statistics
 */
struct ConstantSignalFitnessVals
{
    std::vector<bool> up_states;
    double bursts_seen;
    double norm_max_bursts_diff;
};

/*
 * Function to preprocess a signal and compute statistics
 * Returns ConstantSignalFitnessVals containing up_states, bursts_seen, and norm_max_bursts_diff
 */
ConstantSignalFitnessVals calc_const_signal_vals(const std::vector<double> &signal)
{
    ConstantSignalFitnessVals result;

    double range = HindmarshRose::MAX - HindmarshRose::MIN;
    double th_on = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range + HindmarshRose::MIN;
    double th_up = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range + HindmarshRose::MIN;

    result.up_states.reserve(signal.size());
    bool up = (signal[0] > th_up);
    double bursts_seen = 0;

    for (double val : signal)
    {
        if (!up && val > th_up)
        {
            up = true;
        }
        else if (up && val < th_on)
        {
            up = false;
            bursts_seen++;
        }
        result.up_states.push_back(up);
    }

    result.bursts_seen = bursts_seen;
    result.norm_max_bursts_diff = bursts_seen / 2.0;

    return result;
}

/*
 * Calculate fitness based on bursts activity (XNOR for phase, XOR for antiphase)
 * Takes precomputed stats for signal1, raw signal2, and a bool search_phase (true for phase, false for antiphase)
 * Returns the computed fitness score
 */
double fitness_from_signals(const ConstantSignalFitnessVals &stats1, const std::vector<double> &signal2, bool search_phase)
{
    size_t signal_size = signal2.size();

    double min2 = GeneralConstants::DOUBLE_MAX;
    double max2 = GeneralConstants::DOUBLE_MIN;

    for (double val2 : signal2)
    {
        if (val2 < min2)
            min2 = val2;
        if (val2 > max2)
            max2 = val2;
    }

    // Calculate dynamic thresholds for signal2
    double range2 = max2 - min2;
    double th_on2 = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range2 + min2;
    double th_up2 = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range2 + min2;

    // Compute minmax_score
    double minmax_diff = std::abs(HindmarshRose::MAX - max2) + std::abs(HindmarshRose::MIN - min2);
    double minmax_score = inverse_normalization(minmax_diff, 0.0, FitnessConstants::NORM_MAX_MINMAX_DIFF);

    // State machine for signal2
    bool up2 = (signal2[0] > th_up2);
    double bursts_seen_2 = 0;
    double phase_score = 0.0;

    for (size_t i = 0; i < signal_size; i++)
    {
        double val2 = signal2[i];
        bool up1 = stats1.up_states[i];

        // State machine for signal2
        if (!up2 && val2 > th_up2)
        {
            up2 = true;
        }
        else if (up2 && val2 < th_on2)
        {
            up2 = false;
            bursts_seen_2++;
        }

        // Fitness logic: calculate phase fraction (XNOR)
        if (up1 == up2)
        {
            phase_score += 1.0;
        }
    }

    phase_score /= signal_size;

    // Adjust for antiphase mode
    if (!search_phase)
    {
        phase_score = 1.0 - phase_score;
    }

    // Compute bursts_score using precomputed for signal1
    double bursts_diff_score = inverse_normalization(std::abs(stats1.bursts_seen - bursts_seen_2), 0.0, stats1.norm_max_bursts_diff);

    // Compute final weighted score
    double final_score = (FitnessConfig::BURSTS_DIFF_WEIGHT * bursts_diff_score) + (FitnessConfig::PHASE_WEIGHT * phase_score) + (FitnessConfig::MINMAX_WEIGHT * minmax_score);

    return final_score;
}

/*
 * Template function to calculate fitnesses for multiple parameter sets
 * Simulates the neural model with given parameters and computes fitness against precomputed stats
 * Parameters: synapsis, neurons, params_individuals, scaled_result, initial_state, stats1, search_phase, buffers, reset_model_func, get_v_model_func, set_v_csv_func
 */
template <typename Integrator, typename NeuronType, typename StateType, size_t N>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &csv_neur,
                    NeuronType &model_neur,
                    const std::array<ChemicalSynapsisParams, N> &params_individuals,
                    const ScaledSignalResult &scaled_result,
                    const StateType &initial_state,
                    const ConstantSignalFitnessVals &stats1,
                    bool search_phase,
                    std::vector<double> &model_signal_buffer,
                    std::array<double, N> &fitnesses_buffer,
                    void (*reset_model_func)(NeuronType &, const StateType &),
                    double (*get_v_model_func)(const NeuronType &),
                    void (*set_v_csv_func)(NeuronType &, double))
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    size_t total_size = scaled_result.signal.size() + scaled_result.interpolated_points.size();

    for (size_t i = 0; i < N; i++)
    {
        const ChemicalSynapsisParams &params = params_individuals[i];

        // Set synapsis parameters
        synapsis.set(ChemicalSynapsisType::gfast, params.gfast);
        synapsis.set(ChemicalSynapsisType::Esyn, params.Esyn);
        synapsis.set(ChemicalSynapsisType::sfast, params.sfast);
        synapsis.set(ChemicalSynapsisType::Vfast, params.Vfast);
        synapsis.set(ChemicalSynapsisType::Vslow, params.Vslow);
        synapsis.set(ChemicalSynapsisType::gslow, params.gslow);
        synapsis.set(ChemicalSynapsisType::k1, params.k1);
        synapsis.set(ChemicalSynapsisType::k2, params.k2);
        synapsis.set(ChemicalSynapsisType::sslow, params.sslow);

        // Reset synapsis
        synapsis.set(ChemicalSynapsisType::mslow, FitnessConstants::M_SLOW_INITIAL_VALUE);

        // Reset neuron state using provided function
        reset_model_func(model_neur, initial_state);

        // Simulate and collect neur signal
        size_t signal_counter = 0;
        size_t interpolated_signal_counter = 0;
        for (size_t j = 0; j < total_size; j++)
        {
            double model_neur_val = get_v_model_func(model_neur);

            double csv_neur_val_to_use;
            if ((j % scaled_result.points_factor) == 0)
            {
                csv_neur_val_to_use = scaled_result.signal[signal_counter];
                model_signal_buffer[signal_counter] = model_neur_val;
                signal_counter++;
            }
            else
            {
                csv_neur_val_to_use = scaled_result.interpolated_points[interpolated_signal_counter];
                interpolated_signal_counter++;
            }

            set_v_csv_func(csv_neur, csv_neur_val_to_use);
            synapsis.step(scaled_result.dt, csv_neur_val_to_use, model_neur_val);
            model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
            model_neur.step(scaled_result.dt);
        }

        // Compute fitness
        fitnesses_buffer[i] = fitness_from_signals(stats1, model_signal_buffer, search_phase);
    }
}

/*
 * Reset the state of a Hindmarsh-Rose neuron using the provided StateType
 */
template <typename NeuronType>
void reset_hindmarsh_rose_state(NeuronType &neuron, const HindmarshRoseState &state)
{
    neuron.set(NeuronType::x, state.x);
    neuron.set(NeuronType::y, state.y);
    neuron.set(NeuronType::z, state.z);
    neuron.reset_synaptic_input();
}

/*
 * Get the voltage (membrane potential) from a Hindmarsh-Rose neuron
 */
template <typename NeuronType>
double get_hindmarsh_rose_voltage(const NeuronType &neuron)
{
    return neuron.get(NeuronType::x);
}

/*
 * Set the voltage (membrane potential) of a Hindmarsh-Rose neuron
 */
template <typename NeuronType>
void set_hindmarsh_rose_voltage(NeuronType &neuron, double voltage)
{
    neuron.set(NeuronType::x, voltage);
}