#ifndef SCALING_H
#define SCALING_H

#include <cstddef>
#include <optional>
#include <string>
#include <kfr/all.hpp>
using namespace kfr;

// Numeric integrators supported by the offline pipeline.
enum NumericIntegrator
{
    RK4
};

// Neuron models supported by the offline pipeline.
enum NeuronModel
{
    HINDMARSH_ROSE
};

struct ScaledSigResult
{
    // Input trace rescaled to the model voltage domain.
    univector<double> sig;
    // Intermediate points linearly interpolated between consecutive input samples.
    univector<double> interpolated_points;
    // Number of model integration points per input sample.
    size_t points_factor;
    // Integration step selected for the model simulation (ms).
    double dt;
};

// Read, scale, and resample a CSV voltage signal to the selected model/integrator.
std::optional<ScaledSigResult> scale_sig(
    const std::string &csv_path,
    size_t column_idx,
    double csv_step,
    double start_time,
    double use_time,
    double observation_time,
    NumericIntegrator integrator,
    NeuronModel model,
    bool check_drift);

#endif
