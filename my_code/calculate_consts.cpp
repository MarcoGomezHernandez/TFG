/*
 * calculate_consts.cpp
 * Calcula constantes de normalización (pts, min, max) para modelos neuronales (NeuN).
 */

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <variant>

// NeuN Headers
#include <DifferentialNeuronWrapper.h>
#include <HindmarshRoseModel.h>
#include <SystemWrapper.h>
#include <RungeKutta4.h>
#include "scaling.h"
#include "scaling_utils.h"

// Definición de dts (Input)
static constexpr std::array<double, 144> dts = {
    0.000500, 0.000600, 0.000700, 0.000800, 0.000900, 0.001000, 0.001100, 0.001200,
    0.001300, 0.001400, 0.001500, 0.001600, 0.001800, 0.002000, 0.002200, 0.002500,
    0.002800, 0.002900, 0.003000, 0.003100, 0.003200, 0.003300, 0.003400, 0.003500,
    0.003600, 0.003700, 0.003800, 0.003900, 0.004000, 0.004100, 0.004200, 0.004300,
    0.004400, 0.004500, 0.004600, 0.004700, 0.004800, 0.004900, 0.005000, 0.005100,
    0.005200, 0.005400, 0.005600, 0.005800, 0.006000, 0.006200, 0.006400, 0.006600,
    0.006800, 0.007000, 0.007200, 0.007400, 0.007700, 0.008000, 0.008300, 0.008600,
    0.008900, 0.009200, 0.009600, 0.010000, 0.010400, 0.010900, 0.011400, 0.011900,
    0.012500, 0.013100, 0.013800, 0.014600, 0.015400, 0.016300, 0.017300, 0.018500,
    0.019900, 0.021500, 0.023300, 0.025500, 0.028100, 0.028400, 0.028700, 0.029000,
    0.029400, 0.029800, 0.030200, 0.030600, 0.031000, 0.031400, 0.031800, 0.032200,
    0.032600, 0.033000, 0.033400, 0.033900, 0.034400, 0.034900, 0.035400, 0.035900,
    0.036400, 0.036900, 0.037400, 0.038000, 0.038600, 0.039200, 0.039800, 0.040400,
    0.041000, 0.041700, 0.042400, 0.043100, 0.043800, 0.044500, 0.045300, 0.046100,
    0.046900, 0.047700, 0.048600, 0.049500, 0.050400, 0.051400, 0.052400, 0.053400,
    0.054500, 0.055600, 0.056800, 0.058000, 0.059300, 0.060600, 0.062000, 0.063400,
    0.064900, 0.066500, 0.068200, 0.069900, 0.071700, 0.073600, 0.075600, 0.077700,
    0.079900, 0.082300, 0.084800, 0.087500, 0.090300, 0.093300, 0.096500, 0.100000};

static constexpr size_t dts_size = dts.size();

// Add struct for neuron state
struct HindmarshRoseState
{
    double x;
    double y;
    double z;
};

// Parámetros y estado inicial del modelo Hindmarsh-Rose
struct HindmarshRoseParams
{
    double e;
    double mu;
    double S;
    double a;
    double b;
    double c;
    double d;
    double xr;
    double vh;
};

// Estructura de retorno
struct CalcResult
{
    std::array<double, dts_size> pts;
    double min;
    double max;
    std::vector<double> invalid_dts;
};

namespace ConstCalculatorConstants
{
    static constexpr size_t BURSTS_TO_AVERAGE = 20;
}

/**
 * Función principal de cálculo.
 * N: Tamaño del array dts
 */
template <typename NeuronType, typename NeuronParamsType, typename StateType, size_t N>
CalcResult calculate_metrics(
    NeuronModel model,
    const NeuronParamsType &config,
    const StateType &initial_state,
    const std::array<double, N> &dts,
    double observation_time,
    double minmax_dt,
    double stabilization_time)
{
    if observation_time
        <= 0 || minmax_dt <= 0 || stabilization_time < 0
        {
            throw std::runtime_error("observation_time and minmax_dt must be positive, stabilization_time non-negative");
        }

    NeuronType::ConstructorArgs args;
    if (model == HINDMARSH_ROSE)
    {
        args.params[NeuronType::e] = config.e;
        args.params[NeuronType::mu] = config.mu;
        args.params[NeuronType::S] = config.S;
        args.params[NeuronType::a] = config.a;
        args.params[NeuronType::b] = config.b;
        args.params[NeuronType::c] = config.c;
        args.params[NeuronType::d] = config.d;
        args.params[NeuronType::xr] = config.xr;
        args.params[NeuronType::vh] = config.vh;
    }
    else
    {
        throw std::runtime_error("Modelo no implementado.");
    }

    NeuronType neuron(args);

    CalcResult result;

    // Set to initial state
    if (model == HINDMARSH_ROSE)
    {
        neuron.set(NeuronType::x, initial_state.x);
        neuron.set(NeuronType::y, initial_state.y);
        neuron.set(NeuronType::z, initial_state.z);
    }

    // Stabilization with minmax_dt
    size_t stabilization_steps = static_cast<size_t>(stabilization_time / minmax_dt);
    for (size_t i = 0; i < stabilization_steps; i++)
        neuron.step(minmax_dt);

    double min_abs = SignalConstants::DOUBLE_MAX;
    double max_abs = SignalConstants::DOUBLE_MIN;
    size_t obs_steps = static_cast<size_t>(observation_time / minmax_dt);
    for (size_t step = 0; step < obs_steps; step++)
    {
        neuron.step(minmax_dt);

        double val;
        if (model == HINDMARSH_ROSE)
        {
            val = neuron.get(NeuronType::x);
        }

        if (val > max_abs)
            max_abs = val;
        if (val < min_abs)
            min_abs = val;
    }
    result.min = min_abs;
    result.max = max_abs;

    // Compute relative thresholds
    double range = max_abs - min_abs;
    double th_on = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range + min_abs;
    double th_up = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range + min_abs;

    // Iterar sobre cada dt
    for (size_t i = 0; i < N; i++)
    {
        double dt = dts[i];

        // Set to initial state
        if (model == HINDMARSH_ROSE)
        {
            neuron.set(NeuronType::x, initial_state.x);
            neuron.set(NeuronType::y, initial_state.y);
            neuron.set(NeuronType::z, initial_state.z);
        }

        // Estabilización inicial (transitorio) with current dt
        size_t stabilization_steps = static_cast<size_t>(stabilization_time / dts[i]);
        for (int j = 0; j < stabilization_steps; j++)
            neuron.step(dt);

        // Second pass: detect bursts using relative thresholds
        double total_steps = 0;
        size_t bursts_seen = 0;
        bool up = false;
        size_t steps_in_current_burst = 0;

        // Simulación
        size_t step = 0;
        size_t max_steps = static_cast<size_t>(observation_time / dt);

        while (bursts_seen < ConstCalculatorConstants::BURSTS_TO_AVERAGE && step < max_steps)
        {
            neuron.step(dt);
            step++;

            double val;
            if (model == HINDMARSH_ROSE)
            {
                val = neuron.get(NeuronType::x);
            }

            if (!up && val > th_up)
            {
                up = true;
                steps_in_current_burst = 0;
            }
            else if (up && val < th_on)
            {
                up = false;
                bursts_seen++;
                total_steps += steps_in_current_burst;
            }

            if (up)
            {
                steps_in_current_burst++;
            }
        }

        if (bursts_seen == 0)
        {
            result.pts[i] = SignalConstants::DOUBLE_MAX; // No bursts detected, set to max as sentinel
            result.invalid_dts.push_back(dts[i]);
        }
        else
        {
            result.pts[i] = total_steps / static_cast<double>(bursts_seen);
        }
    }

    return result;
}

int main()
{
    // 2. Configuración
    HindmarshRoseParams hr_params = {
        1.0,   // a
        3.0,   // b
        1.0,   // c
        5.0,   // d
        0.006, // r
        4.0,   // s
        -1.6,  // x_rest
        3.0    // I (genera bursting)
    };

    HindmarshRoseState initial_state = {
        -1.6,  // x
        -10.0, // y
        0.0    // z
    };

    // 3. Ejecución
    double observation_time = 1000.0;   // Dimensionless observation time
    double minmax_dt = dts[0];          // Use first dt for min/max calculation
    double stabilization_time = 1000.0; // Stabilization time
    std::cout << "Calculando constantes... Espere.\n";
    auto result = calculate_metrics<DifferentialNeuronWrapper<SystemWrapper<HindmarshRoseModel<double>>, RungeKutta4>, HindmarshRoseParams, HindmarshRoseState>(HINDMARSH_ROSE, hr_params, initial_state, dts, observation_time, minmax_dt, stabilization_time);

    // 4. Salida Formateada
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\n// Copiar y pegar en scaling.cpp / struct correspondiente:\n\n";
    std::cout << "static constexpr double min = " << result.min << ";\n";
    std::cout << "static constexpr double max = " << result.max << ";\n";

    std::cout << "static constexpr std::array<double, " << dts_size << "> dts = {";
    for (size_t i = 0; i < dts_size; ++i)
    {
        if (i % 8 == 0)
            std::cout << "\n    ";
        std::cout << dts[i];
        if (i < dts_size - 1)
            std::cout << ", ";
    }
    std::cout << "};\n";

    std::cout << "static constexpr std::array<double, " << dts_size << "> pts = {";
    for (size_t i = 0; i < dts_size; ++i)
    {
        if (i % 8 == 0)
            std::cout << "\n    ";
        std::cout << result.pts[i];
        if (i < dts_size - 1)
            std::cout << ", ";
    }
    std::cout << "};\n";

    if (!result.invalid_dts.empty())
    {
        std::cout << "Invalid dts (no bursts detected): ";
        for (auto dt : result.invalid_dts)
        {
            std::cout << dt << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}