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

/*
 * Calculate fitness based on exclusive burst activity (XOR logic)
 * Returns a score based on mismatched burst states, or 0 if signals are invalid
 */
double antiphase_fitness(const std::vector<double> &signal1, const std::vector<double> &signal2)
{
    size_t signal_size = signal1.size();

    double min1 = SignalConstants::DOUBLE_MAX;
    double max1 = SignalConstants::DOUBLE_MIN;
    double min2 = SignalConstants::DOUBLE_MAX;
    double max2 = SignalConstants::DOUBLE_MIN;

    for (size_t i = 0; i < signal_size; i++)
    {
        if (signal1[i] < min1)
            min1 = signal1[i];
        if (signal1[i] > max1)
            max1 = signal1[i];

        if (signal2[i] < min2)
            min2 = signal2[i];
        if (signal2[i] > max2)
            max2 = signal2[i];
    }

    // Calculate dynamic thresholds
    double range1 = max1 - min1;
    double th_on1 = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range1 + min1;
    double th_up1 = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range1 + min1;

    double range2 = max2 - min2;
    double th_on2 = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range2 + min2;
    double th_up2 = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range2 + min2;

    // 2. Second Pass: State machine for burst detection and fitness accumulation
    bool up1 = (signal1[0] > th_up1);
    bool up2 = (signal2[0] > th_up2);

    int bursts_seen1 = 0;
    int bursts_seen2 = 0;

    double antiphase_score = 0.0;

    for (size_t i = 0; i < signal_size; i++)
    {
        double val1 = signal1[i];
        double val2 = signal2[i];

        // --- State Machine Signal 1 ---
        if (!up1 && val1 > th_up1)
        {
            up1 = true;
        }
        else if (up1 && val1 < th_on1)
        {
            up1 = false;
            bursts_seen1++;
        }

        // --- State Machine Signal 2 ---
        if (!up2 && val2 > th_up2)
        {
            up2 = true;
        }
        else if (up2 && val2 < th_on2)
        {
            up2 = false;
            bursts_seen2++;
        }

        // --- Fitness Logic ---
        // Increment score if one signal is in burst and the other is not (XOR)
        if (up1 != up2)
        {
            antiphase_score += 1.0;
        }
    }

    antiphase_score /= signal_size; // Normalize to [0,1]

    // Compute bursts_score
    double burst_min_score_1 = min_0_no_max_desc_normalization(bursts_seen1, FitnessConfig::MINIMUM_BURSTS_WITH_MAX_SCORE);
    double burst_min_score_2 = min_0_no_max_desc_normalization(bursts_seen2, FitnessConfig::MINIMUM_BURSTS_WITH_MAX_SCORE);
    double burst_diff_score = hard_sigmoid(std::abs(bursts_seen1 - bursts_seen2), FitnessConfig::HALF_SCORE_BURST_DIFFERENCE);
    double bursts_score = (FitnessConfig::BURST_MIN_WEIGHT * burst_min_score_1) +
                          (FitnessConfig::BURST_MIN_WEIGHT * burst_min_score_2) +
                          (FitnessConfig::BURST_DIFF_WEIGHT * burst_diff_score);

    // Compute final weighted score
    double final_score = (FitnessConfig::BURSTS_WEIGHT * bursts_score) + (FitnessConfig::ANTIPHASE_WEIGHT * antiphase_score);

    return final_score;
}