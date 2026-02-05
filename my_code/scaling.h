#ifndef SCALING_H
#define SCALING_H

#include <vector>
#include <string>

// Constants
constexpr int RK4 = 0;
constexpr double HINDMARSH_ROSE_MIN = -1.608734;
constexpr double HINDMARSH_ROSE_MAX = 1.797032;

// Structure to hold the result
struct ScaledSignalResult {
    std::vector<double> scaled_signal;
    double dt;
    bool success;
};

// Main external function
ScaledSignalResult scale_signal(
    const std::string& csv_path,
    size_t column_index,
    double csv_step,
    double start_time,
    double use_time,
    double observation_time,
    double freq,
    int integrator,
    bool check_drift,
    double sec_per_burst = -1.0
);

#endif // SCALING_H
