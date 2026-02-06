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

// Variant for config types
using NeuronParams = std::variant<HindmarshRoseParams>;

static constexpr size_t dts_size = dts.size();

// Estructura de retorno
struct CalcResult
{
    std::array<double, dts_size> pts;
    double global_min;
    double global_max;
};

/**
 * Función principal de cálculo.
 * N: Tamaño del array dts
 */
template <size_t N>
CalcResult calculate_metrics(
    NumericIntegrator integrator,
    NeuronModel model,
    const NeuronParams &config,
    const std::array<double, N> &dts,
    size_t periods_to_average = 10)
{
    CalcResult<N> result;
    result.global_min = std::numeric_limits<double>::max();
    result.global_max = std::numeric_limits<double>::lowest();

    // Lógica para Hindmarsh-Rose con RK4
    if (model == HINDMARSH_ROSE)
    {

        if (integrator == RK4)
        {
            using Integrator = RungeKutta4;
        }
        else
        {
            throw std::runtime_error("Integrador no implementado para Hindmarsh-Rose.");
        }
        using HR_Type = DifferentialNeuronWrapper<SystemWrapper<HindmarshRoseModel<double>>, Integrator>;

        // Mapeo de parámetros
        // HR Param order: a, b, c, d, r, s, x_rest, I
        HR_Type::ConstructorArgs args;
        if (std::holds_alternative<HindmarshRoseParams>(config))
        {
            auto &p = std::get<HindmarshRoseParams>(config);
            args.params[0] = p.a;
            args.params[1] = p.b;
            args.params[2] = p.c;
            args.params[3] = p.d;
            args.params[4] = p.r;
            args.params[5] = p.s;
            args.params[6] = p.x_rest;
            args.params[7] = p.I;
        }

        // Iterar sobre cada dt
        for (size_t i = 0; i < N; ++i)
        {
            double dt = dts[i];

            // Instanciar modelo
            HR_Type neuron(args);
            if (std::holds_alternative<HindmarshRoseParams>(config))
            {
                auto &p = std::get<HindmarshRoseParams>(config);
                neuron.set(HR_Type::x, p.x);
                neuron.set(HR_Type::y, p.y);
                neuron.set(HR_Type::z, p.z);
            }

            // Variables para detección de periodo
            int burst_count = 0;
            bool in_burst = false;
            double threshold = 1.0; // Umbral empírico para detección de disparo en HR

            double time_start = 0.0;
            double time_end = 0.0;
            double current_time = 0.0;

            // Simulación
            // Límite de seguridad: 1 millón de pasos o suficientes periodos
            size_t step = 0;
            size_t max_steps_safety = 20000000;

            // Estabilización inicial (transitorio)
            for (int k = 0; k < 10000; ++k)
                neuron.step(dt);

            while (burst_count <= periods_to_average && step < max_steps_safety)
            {
                double x_val = neuron.get(HR_Type::x);

                // Actualizar min/max globales
                if (x_val < result.global_min)
                    result.global_min = x_val;
                if (x_val > result.global_max)
                    result.global_max = x_val;

                // Detección de flanco de subida (inicio de burst)
                if (x_val > threshold && !in_burst)
                {
                    in_burst = true;
                    if (burst_count == 0)
                    {
                        time_start = current_time;
                    }
                    if (burst_count == periods_to_average)
                    {
                        time_end = current_time;
                    }
                    burst_count++;
                }
                else if (x_val < 0.0 && in_burst)
                { // Hysteresis reset
                    in_burst = false;
                }

                neuron.step(dt);
                current_time += dt;
                step++;
            }

            if (burst_count > periods_to_average)
            {
                double total_time = time_end - time_start;
                double avg_period = total_time / static_cast<double>(periods_to_average);
                result.pts[i] = avg_period / dt;
            }
            else
            {
                std::cerr << "Warning: No se detectaron suficientes periodos para dt=" << dt << "\n";
                result.pts[i] = 0.0;
            }
        }
    }
    else
    {
        throw std::runtime_error("Modelo no implementado.");
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
        3.0,   // I (genera bursting)
        -1.6,  // x
        -10.0, // y
        0.0    // z
    };

    // 3. Ejecución
    std::cout << "Calculando constantes... Espere.\n";
    auto result = calculate_metrics(RK4, HINDMARSH_ROSE, hr_params, dts, 20);

    // 4. Salida Formateada
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\n// Copiar y pegar en scaling.cpp / struct correspondiente:\n\n";
    std::cout << "static constexpr double min = " << result.global_min << ";\n";
    std::cout << "static constexpr double max = " << result.global_max << ";\n";

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

    return 0;
}