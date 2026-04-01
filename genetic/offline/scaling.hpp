#ifndef SCALING_H
#define SCALING_H

#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <iostream>
#include <kfr/all.hpp>
using namespace kfr;
#include "utils.hpp"

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
    bool success;
};

ScaledSigResult scale_sig(
    const std::string &csv_path,
    size_t column_i,
    double csv_step,
    double start_time,
    double use_time,
    double observation_time,
    NumericIntegrator integrator,
    NeuronModel model,
    bool check_drift);

#endif
