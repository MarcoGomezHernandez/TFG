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
    size_t signal_size = signal1.size();

    // Validation: Sizes must match
    if (signal_size != signal2.size() || signal_size == 0)
    {
        return 0.0;
    }

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
            antiphase_score += 1.0;
        }
    }

    // 3. Validation: Check minimum burst requirements
    if (bursts_seen1 < FitnessConfig::MIN_REQUIRED_BURSTS ||
        bursts_seen2 < FitnessConfig::MIN_REQUIRED_BURSTS)
    {
        return 0.0;
    }

    return antiphase_score;
}