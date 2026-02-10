#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "signal_utils.h"

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

/*
 * Compute an inverse normalization that maps a value to a score between 0 and 1, with score=1 at val=min_val and score=0 at val=max_val, clamped to [0,1].
 */
double inverse_normalization(double val, double min_val, double max_val)
{
    return (max_val - val) / (max_val - min_val);
}

// struct to hold precomputed signal statistics
struct ConstantSignalFitnessVals
{
    double min;
    double max;
    std::vector<bool> up_states;
    double bursts_seen;
    double norm_max_bursts_diff;
    double norm_max_minmax_diff;
};

// function to preprocess a signal and compute statistics
ConstantSignalFitnessVals calc_const_signal_vals(const std::vector<double> &signal, double min, double max)
{
    ConstantSignalFitnessVals result;

    result.min = min;
    result.max = max;

    double range = max - min;
    double th_on = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range + min;
    double th_up = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range + min;

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
    result.norm_max_minmax_diff = (HindmarshRose::MAX - HindmarshRose::MIN) * 2.0;

    return result;
}

/*
 * Calculate fitness based on bursts activity (XNOR for phase, XOR for antiphase)
 * Now takes precomputed stats for signal1, raw signal2, and a bool search_phase (true for phase, false for antiphase)
 */
double fitness_from_signals(const ConstantSignalFitnessVals &stats1, const std::vector<double> &signal2, bool search_phase)
{
    size_t signal_size = signal2.size();

    double min2 = SignalConstants::DOUBLE_MAX;
    double max2 = SignalConstants::DOUBLE_MIN;

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
    double minmax_diff = std::abs(stats1.max - max2) + std::abs(stats1.min - min2);
    double minmax_score = inverse_normalization(minmax_diff, 0.0, stats1.norm_max_minmax_diff);

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