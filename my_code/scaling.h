#ifndef SCALING_H
#define SCALING_H

#include <vector>
#include <string>

// Integrator types
enum Integrator {
    RK4
};

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
    Integrator integrator,
    bool check_drift
);

#endif // SCALING_H
