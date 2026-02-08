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
    // Minimum bursts required to consider the signal valid
    static constexpr int MIN_REQUIRED_BURSTS = 3;
}

/*
 * Calculate fitness based on exclusive burst activity (XOR logic)
 * Returns a score based on mismatched burst states, or 0 if signals are invalid
 */
double fitness(const std::vector<double> &signal1, const std::vector<double> &signal2)
{
    // Validation: Sizes must match
    if (signal1.size() != signal2.size() || signal1.empty())
    {
        return 0.0;
    }

    size_t n_steps = signal1.size();

    // 1. First Pass: Calculate Min/Max for both signals
    auto [min1_it, max1_it] = std::minmax_element(signal1.begin(), signal1.end());
    auto [min2_it, max2_it] = std::minmax_element(signal2.begin(), signal2.end());

    double min1 = *min1_it;
    double max1 = *max1_it;
    double min2 = *min2_it;
    double max2 = *max2_it;

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

    double fitness_score = 0.0;

    for (size_t i = 0; i < n_steps; i++)
    {
        double val1 = signal1[i];
        double val2 = signal2[i];

        // --- State Machine Signal 1 ---
        if (!up1 && val1 > th_up1)
        {
            up1 = true;
            bursts_seen1++;
        }
        else if (up1 && val1 < th_on1)
        {
            up1 = false;
        }

        // --- State Machine Signal 2 ---
        if (!up2 && val2 > th_up2)
        {
            up2 = true;
            bursts_seen2++;
        }
        else if (up2 && val2 < th_on2)
        {
            up2 = false;
        }

        // --- Fitness Logic ---
        // Increment score if one signal is in burst and the other is not (XOR)
        if (up1 != up2)
        {
            fitness_score += 1.0;
        }
    }

    // 3. Validation: Check minimum burst requirements
    if (bursts_seen1 < FitnessConfig::MIN_REQUIRED_BURSTS ||
        bursts_seen2 < FitnessConfig::MIN_REQUIRED_BURSTS)
    {
        return 0.0;
    }

    return fitness_score;
}