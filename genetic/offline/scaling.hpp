#ifndef SCALING_H
#define SCALING_H

#include <cstddef>
#include <optional>
#include <string>
#include <kfr/all.hpp>
using namespace kfr;

enum NumericIntegrator
{
    RK4
};

enum NeuronModel
{
    HINDMARSH_ROSE
};

struct ScaledSigResult
{
    univector<double> sig;
    univector<double> interpolated_points;
    size_t points_factor;
    double dt;
};

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
