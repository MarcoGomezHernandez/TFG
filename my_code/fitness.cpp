#include "fitness.hpp"

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
 * Function to preprocess a signal and compute statistics
 * Returns ConstantSignalFitnessVals containing up_states, bursts_seen, and norm_max_bursts_diff
 */
ConstantSignalFitnessVals calc_const_signal_vals(const std::vector<double> &signal, double min_val, double max_val)
{
    ConstantSignalFitnessVals result;

    const double range = max_val - min_val;
    const double th_on = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range + min_val;
    const double th_up = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range + min_val;

    std::vector<bool> &up_states = result.up_states;
    up_states.reserve(signal.size());
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
        up_states.push_back(up);
    }

    result.bursts_seen = bursts_seen;
    result.norm_max_bursts_diff = bursts_seen / 2.0;
    result.min_val = min_val;
    result.max_val = max_val;

    return result;
}

/*
 * Calculate fitness based on bursts activity (XNOR for phase, XOR for antiphase)
 * Takes precomputed stats for signal1, raw signal2, and a bool search_phase (true for phase, false for antiphase)
 * Returns the computed fitness score
 */
double fitness_from_signals(const ConstantSignalFitnessVals &stats1, const std::vector<double> &signal2, bool search_phase)
{
    const size_t signal_size = signal2.size();

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
    const double range2 = max2 - min2;
    const double th_on2 = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range2 + min2;
    const double th_up2 = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range2 + min2;

    // Compute minmax_score
    const double minmax_diff = std::abs(stats1.max_val - max2) + std::abs(stats1.min_val - min2);
    const double minmax_score = inverse_normalization(minmax_diff, 0.0, FitnessConstants::NORM_MAX_MINMAX_DIFF);

    // State machine for signal2
    bool up2 = (signal2[0] > th_up2);
    double bursts_seen_2 = 0;
    double phase_score = 0.0;

    const std::vector<bool> &up_states = stats1.up_states;
    for (size_t i = 0; i < signal_size; i++)
    {
        double val2 = signal2[i];
        bool up1 = up_states[i];

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
    const double bursts_diff_score = inverse_normalization(std::abs(stats1.bursts_seen - bursts_seen_2), 0.0, stats1.norm_max_bursts_diff);

    // Compute final weighted score
    const double final_score = (FitnessConfig::BURSTS_DIFF_WEIGHT * bursts_diff_score) + (FitnessConfig::PHASE_WEIGHT * phase_score) + (FitnessConfig::MINMAX_WEIGHT * minmax_score);

    return final_score;
}

/*
 * Template function to calculate fitnesses for multiple parameter sets
 * Simulates the neural model with given parameters and computes fitness against precomputed stats
 * Parameters: synapsis, neurons, individuals, scaled_result, stats1, search_phase, buffers, reset_state_neur, get_v_neur
 */
template <typename Integrator, typename NeuronType, size_t N, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    std::array<Individual, N> &individuals,
                    const ScaledSignalResult &scaled_result,
                    const ConstantSignalFitnessVals &stats1,
                    bool search_phase,
                    std::vector<double> &model_signal_buffer,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur)
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    const size_t signal_size = scaled_result.signal.size();
    const size_t points_factor = scaled_result.points_factor;
    const double dt = scaled_result.dt;
    const double *signal_data = scaled_result.signal.data();
    const double *interpolated_signal_data = scaled_result.interpolated_signal.data();

    for (size_t i = 0; i < N; i++)
    {
        const ChemicalSynapsisParams &params = individuals[i].params;

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
        reset_state_neur(model_neur);

        // Simulate and collect neur signal
        size_t interp_signal_counter = 0;
        size_t j = 0;
        for (; j < signal_size; j++)
        {
            synapsis.step(dt, signal_data[j], get_v_neur(model_neur));
            model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
            model_neur.step(dt);
            model_signal_buffer[j] = get_v_neur(model_neur);

            for (size_t k = 1; k < points_factor; k++)
            {
                synapsis.step(dt, interpolated_signal_data[interp_signal_counter], get_v_neur(model_neur));
                model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
                model_neur.step(dt);
                interp_signal_counter++;
            }
        }

        synapsis.step(dt, signal_data[j], get_v_neur(model_neur));
        model_neur.add_synaptic_input(synapsis.get(ChemicalSynapsisType::i));
        model_neur.step(dt);
        model_signal_buffer[j] = get_v_neur(model_neur);

        // Compute fitness and store directly in individual
        individuals[i].fitness = fitness_from_signals(stats1, model_signal_buffer, search_phase);
    }
}