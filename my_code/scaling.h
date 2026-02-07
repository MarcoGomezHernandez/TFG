#ifndef SCALING_H
#define SCALING_H

#include <vector>
#include <string>

namespace SignalConstants
{
    // Sentinel values
    inline constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();
    inline constexpr double DOUBLE_MIN = std::numeric_limits<double>::lowest();

    // Invalid values
    inline constexpr double INVALID_DT = -1.0;
    inline constexpr double INVALID_PTS = -1.0;
}

// Integrator types
enum NumericIntegrator
{
    RK4
};

// Neuron model types
enum NeuronModel
{
    HINDMARSH_ROSE
};

// Structure to hold the result
struct ScaledSignalResult
{
    std::vector<double> scaled_signal;
    double dt;
    bool success;
};

// Main external function
ScaledSignalResult scale_signal(
    const std::string &csv_path,
    size_t column_index,
    double csv_step,
    double start_time,
    double use_time,
    double observation_time,
    NumericIntegrator integrator,
    NeuronModel model,
    bool check_drift);

#endif // SCALING_H
