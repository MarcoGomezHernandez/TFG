// Shared constants and configuration for signal scaling
#ifndef SCALING_UTILS_H
#define SCALING_UTILS_H

#include <limits>

/*
 * Sentinel and invalid values used throughout signal processing
 */
namespace SignalConstants
{
    // Sentinel values for range tracking
    inline constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();
    inline constexpr double DOUBLE_MIN = std::numeric_limits<double>::lowest();

    // Invalid/uninitialized values
    inline constexpr double INVALID_DT = -1.0;
    inline constexpr double INVALID_PTS = -1.0;
}

/*
 * Public configuration for signal processing
 */
namespace SignalPublicConfig
{
    // Relative thresholds for period detection (10%-90% of signal range)
    inline constexpr double SIGNAL_PERCENTAGE_MIN = 0.10;
    inline constexpr double SIGNAL_PERCENTAGE_MAX = 0.90;
}

#endif // SCALING_UTILS_H
