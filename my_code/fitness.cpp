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
    // Weight for scores
    static constexpr double BURSTS_WEIGHT = 0.5;
    static constexpr double ANTIPHASE_WEIGHT = 0.5;
    static constexpr double MINIMUM_BURSTS_WITH_MAX_SCORE = 20.0; // Minimum number of bursts to achieve maximum score in burst count component. We assume that in the perfect case we will see 20 bursts
    // Sensibility constant for algebraic sigmoid; set to the value that will get half (0.5) of the score.
    static constexpr double HALF_SCORE_BURST_DIFFERENCE = 1.0; // Very low value to heavily penalize any difference in burst count between the two signals
    // Weights for burst score components
    static constexpr double BURST_MIN_WEIGHT = 0.25;
    static constexpr double BURST_DIFF_WEIGHT = 0.5;
}

/*
 * Compute a function that maps a non-negative value to a score between 0 and 1, with a maximum score of 1 at val=0 and approaching 0 as val increases. The sensibility parameter controls how quickly the score decreases as val increases, and must be positive.
 */
double min_0_no_max_desc_normalization(double val, double sensibility)
{
    return 1.0 - (val / (val + sensibility));
}

/*
 * Compute a function that maps a non-negative value to a score between 0 and 1, with a maximum score of 1 at val=min_one_val and the next, a minimum score of 0 at val=0, and a linear score between 0 and 1 for values between 0 and min_one_val. val must be non-negative and min_one_val must non zero.
 */
double hard_sigmoid(double val, double min_one_val)
{
    if (val <= 0)
        return 0.0;
    else if (val >= min_one_val)
        return 1.0;
    else
        return val / min_one_val;
}

// New struct to hold precomputed signal statistics
struct ConstantSignalFitnessVals
{
    double min;
    double max;
    std::vector<bool> up_states;
    double bursts_seen;
    double burst_min_score;
};

// New function to preprocess a signal and compute statistics
ConstantSignalFitnessVals calc_constant_signal_fitness_vals(const std::vector<double> &signal)
{
    ConstantSignalFitnessVals result;

    double min = SignalConstants::DOUBLE_MAX;
    double max = SignalConstants::DOUBLE_MIN;

    for (double val : signal)
    {
        if (val < min)
            min = val;
        if (val > max)
            max = val;
    }

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

    double burst_min_score = min_0_no_max_desc_normalization(bursts_seen, FitnessConfig::MINIMUM_BURSTS_WITH_MAX_SCORE);

    result.bursts_seen = bursts_seen;
    result.burst_min_score = burst_min_score;

    return result;
}

/*
 * Calculate fitness based on exclusive burst activity (XOR logic)
 * Now takes precomputed stats for signal1 and raw signal2
 */
double antiphase_fitness(const ConstantSignalFitnessVals &stats1, const std::vector<double> &signal2)
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

    // State machine for signal2
    bool up2 = (signal2[0] > th_up2);
    double bursts_seen_2 = 0;
    double antiphase_score = 0.0;

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

        // Fitness logic: XOR
        if (up1 != up2)
        {
            antiphase_score += 1.0;
        }
    }

    antiphase_score /= signal_size;

    // Compute bursts_score using precomputed for signal1
    double burst_min_score_2 = min_0_no_max_desc_normalization(bursts_seen_2, FitnessConfig::MINIMUM_BURSTS_WITH_MAX_SCORE);
    double burst_diff_score = hard_sigmoid(std::abs(stats1.bursts_seen - bursts_seen_2), FitnessConfig::HALF_SCORE_BURST_DIFFERENCE);
    double bursts_score = (FitnessConfig::BURST_MIN_WEIGHT * stats1.burst_min_score) +
                          (FitnessConfig::BURST_MIN_WEIGHT * burst_min_score_2) +
                          (FitnessConfig::BURST_DIFF_WEIGHT * burst_diff_score);

    // Compute final weighted score
    double final_score = (FitnessConfig::BURSTS_WEIGHT * bursts_score) + (FitnessConfig::ANTIPHASE_WEIGHT * antiphase_score);

    return final_score;
}