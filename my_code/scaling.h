// Signal scaling interface for neural model integration
#ifndef SCALING_H
#define SCALING_H

#include <vector>
#include <string>

/*
 * Numeric integration methods
 */
enum NumericIntegrator
{
    RK4
};

/*
 * Supported neural models
 */
enum NeuronModel
{
    HINDMARSH_ROSE
};

/*
 * Result of signal scaling operation
 */
struct ScaledSignalResult
{
    std::vector<double> scaled_signal;
    double dt;
    bool success;
};

/*
 * Scale external signal to match neural model dynamics
 */
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
