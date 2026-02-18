
#include <algorithm>
#include <cmath>
#include <iostream>

/*
 * Configuration constants for fitness calculation
 */
namespace FitnessConfig
{
    // Weights for score components (must sum to 1.0)
    static constexpr double BURSTS_DIFF_WEIGHT = 0.1;
    static constexpr double PHASE_WEIGHT = 0.4;
    static constexpr double MINMAX_WEIGHT = 0.1;
    static constexpr double BIOLOGICAL_SIMILARITY_WEIGHT = 0.4;
}

namespace FitnessConstants
{
    static constexpr double NORM_MAX_MINMAX_DIFF = (HindmarshRose::MAX - HindmarshRose::MIN) * 2.0;
    static constexpr double M_SLOW_INITIAL_VALUE = 0.0;
}

/*
 * Compute an inverse normalization that maps a value to a score between 0 and 1, with score=1 at val=min_val and score=0 at val=max_val, clamped to [0,1].
 */
inline double inverse_normalization(double val, double min_val, double max_val)
{
    return (max_val - val) / (max_val - min_val);
}

/*
 * Function to preprocess a signal and compute statistics
 * Returns ConstantSignalFitnessVals containing model_up_states, normalized_signal, bursts_seen, and norm_max_bursts_diff
 * search_phase: true for phase (normalized as-is), false for antiphase (normalized and y-inverted)
 */
ConstantSignalFitnessVals calc_const_signal_fitness_vals(const std::vector<double> &signal, double min_val, double max_val, bool search_phase)
{
    ConstantSignalFitnessVals result;

    const double range = max_val - min_val;
    const double th_on = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range + min_val;
    const double th_up = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range + min_val;

    std::vector<bool> &up_states = result.model_up_states;
    std::vector<double> &normalized_signal = result.normalized_signal;
    const size_t signal_size = signal.size();
    up_states.reserve(signal_size);
    normalized_signal.reserve(signal_size);
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

        // Normalize value to [0, 1] using the model's min/max
        normalized_signal.push_back((val - min_val) / range);
    }

    // For antiphase, invert the normalized signal along the y-axis in a single pass
    if (!search_phase)
    {
        for (double &val : normalized_signal)
            val = 1.0 - val;
    }

    result.bursts_seen = bursts_seen;
    result.norm_max_bursts_diff = bursts_seen / 2.0;
    result.min_val = min_val;
    result.max_val = max_val;

    return result;
}

/*
 * Calculate fitness based on bursts activity (XNOR for phase, XOR for antiphase)
 * Takes precomputed stats for signal1, raw signal2, synapsis_signal, and a bool search_phase (true for phase, false for antiphase)
 * Returns the computed fitness score
 */
double fitness_from_signals(const ConstantSignalFitnessVals &living_stats, const std::vector<double> &model_signal, bool search_phase, const std::vector<double> &synapsis_signal)
{
    const size_t signal_size = model_signal.size();

    double model_min = GeneralConstants::DOUBLE_MAX;
    double model_max = GeneralConstants::DOUBLE_MIN;

    for (double model_val : model_signal)
    {
        if (model_val < model_min)
            model_min = model_val;
        if (model_val > model_max)
            model_max = model_val;
    }

    // Calculate dynamic thresholds for model signal
    const double model_range = model_max - model_min;
    const double model_th_on = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * model_range + model_min;
    const double model_th_up = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * model_range + model_min;

    // Compute minmax_score
    const double minmax_diff = std::abs(living_stats.max_val - model_max) + std::abs(living_stats.min_val - model_min);
    const double minmax_score = inverse_normalization(minmax_diff, 0.0, FitnessConstants::NORM_MAX_MINMAX_DIFF);

    // State machine for model signal
    bool model_up = (model_signal[0] > model_th_up);
    double model_bursts_seen = 0;
    double phase_score = 0.0;

    const std::vector<bool> &up_states = living_stats.model_up_states;
    for (size_t i = 0; i < signal_size; i++)
    {
        double model_val = model_signal[i];
        bool living_up = up_states[i];

        // State machine for model signal
        if (!model_up && model_val > model_th_up)
        {
            model_up = true;
        }
        else if (model_up && model_val < model_th_on)
        {
            model_up = false;
            model_bursts_seen++;
        }

        // Fitness logic: calculate phase fraction (XNOR)
        if (living_up == model_up)
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

    // Compute bursts_score using precomputed for living signal
    const double bursts_diff_score = inverse_normalization(std::abs(living_stats.bursts_seen - model_bursts_seen), 0.0, living_stats.norm_max_bursts_diff);

    // Compute biological_similarity_score: compare normalized reference signal vs normalized synapsis signal
    double syn_min = GeneralConstants::DOUBLE_MAX;
    double syn_max = GeneralConstants::DOUBLE_MIN;
    for (double syn_val : synapsis_signal)
    {
        if (syn_val < syn_min)
            syn_min = syn_val;
        if (syn_val > syn_max)
            syn_max = syn_val;
    }
    const double syn_range = syn_max - syn_min;
    double bio_sim_dist = 0.0;
    const std::vector<double> &living_norm_signal = living_stats.normalized_signal;
    for (size_t i = 0; i < signal_size; i++)
    {
        const double norm_syn_val = (synapsis_signal[i] - syn_min) / syn_range;
        bio_sim_dist += std::abs(living_norm_signal[i] - norm_syn_val);
    }
    // Normalize: max possible accumulated distance is 1.0 * signal_size (both signals in [0,1])
    const double bio_sim_score = 1.0 - (bio_sim_dist / signal_size);

    // Compute final weighted score
    const double final_score = (FitnessConfig::BURSTS_DIFF_WEIGHT * bursts_diff_score) + (FitnessConfig::PHASE_WEIGHT * phase_score) + (FitnessConfig::MINMAX_WEIGHT * minmax_score) + (FitnessConfig::BIOLOGICAL_SIMILARITY_WEIGHT * bio_sim_score);

    if (final_score < 0.0)
    {
        return 0.0;
    }
    return final_score;
}

/*
 * Template function to calculate fitnesses for multiple parameter sets
 * Simulates the neural model with given parameters and computes fitness against precomputed stats
 * Parameters: synapsis, neurons, individuals, scaled_result, stats1, search_phase, buffers, reset_state_neur, get_v_neur, start_index
 */
template <typename Integrator, typename NeuronType, size_t N, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    std::array<Individual, N> &individuals,
                    const ScaledSignalResult &scaled_result,
                    const ConstantSignalFitnessVals &stats1,
                    bool search_phase,
                    std::vector<double> &model_signal_buffer,
                    std::vector<double> &synapsis_signal_buffer,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur,
                    size_t start_index)
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
        for (; j < signal_size; j++)
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
        individuals[i].fitness = fitness_from_signals(stats1, model_signal_buffer, search_phase, synapsis_signal_buffer);
    }
}