#ifndef SCALING_UTILS_H
#define SCALING_UTILS_H

#include <limits>

namespace SignalConstants
{
    // Sentinel values
    inline constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();
    inline constexpr double DOUBLE_MIN = std::numeric_limits<double>::lowest();

    // Invalid values
    inline constexpr double INVALID_DT = -1.0;
    inline constexpr double INVALID_PTS = -1.0;
}

namespace SignalPublicConfig
{
    // Signal range percentages
    static constexpr double SIGNAL_PERCENTAGE_MIN = 0.10;
    static constexpr double SIGNAL_PERCENTAGE_MAX = 0.90;
}

namespace

#endif // SCALING_UTILS_H
