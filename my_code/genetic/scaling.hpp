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
#include "utils.hpp"

enum NumericIntegrator
{
    RK4
};

enum NeuronModel
{
    HINDMARSH_ROSE
};

struct ScaledSignalResult
{
    kfr::univector<double> signal;
    kfr::univector<double> interpolated_points;
    size_t points_factor;
    double dt;
    double pts_burst_real;
    bool success;
};

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

#endif
