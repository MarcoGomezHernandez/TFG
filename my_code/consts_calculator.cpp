// Calculate normalization constants (pts, min, max) for neural models
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
#include "signal_utils.h"

/*
 * State variables for Hindmarsh-Rose model
 */
struct HindmarshRoseState
{
    double x; // Membrane potential
    double y; // Recovery variable
    double z; // Slow adaptation current
};

/*
 * Parameters for Hindmarsh-Rose model configuration
 */
struct HindmarshRoseParams
{
    double e;  // Time scale parameter
    double mu; // Slow dynamics parameter
    double S;  // External stimulus
    double a;  // Cubic nonlinearity coefficient
    double b;  // Quadratic coefficient
    double c;  // Recovery variable coefficient
    double d;  // Recovery variable coefficient
    double xr; // Rest potential
    double vh; // Threshold parameter
};

/*
 * Algorithm constants for burst detection
 */
namespace ConstCalculatorConstants
{
    // Number of bursts to average for accurate pts calculation
    static constexpr int BURSTS_TO_AVERAGE = 20;
}

/*
 * Configuration values for constant calculation
 * Contains dt array, model parameters, and simulation settings
 */
namespace ConstCalculatorConfig
{
    // Array of time steps to test
    static constexpr std::array<double, 144> DTS = {
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

    // Hindmarsh-Rose model parameters
    static constexpr HindmarshRoseParams HR_PARAMS = {
        3.281,  // e
        0.0021, // mu
        1.0,    // S
        1.0,    // a
        3.0,    // b
        1.0,    // c
        5.0,    // d
        -1.6,   // xr
        0.1     // vh
    };

    // Initial state for simulation
    static constexpr HindmarshRoseState HR_INITIAL_STATE = {
        -0.712841, // x
        -1.93688,  // y
        3.16568    // z
    };

    // Simulation time parameters
    static constexpr double OBSERVATION_TIME = 2000.0;   // Time to observe bursts
    static constexpr double MINMAX_DT = DTS[0];          // Finest dt for min/max calculation
    static constexpr double STABILIZATION_TIME = 2000.0; // Transient settling time
}

static constexpr size_t DTS_SIZE = ConstCalculatorConfig::DTS.size();

/*
 * Result structure for min and max values
 */
struct MinMaxResult
{
    double min; // Minimum model output value
    double max; // Maximum model output value
};

/*
 * Result structure for points per burst and invalid dts
 */
struct PtsResult
{
    std::array<double, DTS_SIZE> pts; // Points per burst for each dt
    std::vector<double> invalid_dts;  // dts where no bursts were detected
};

/*
 * Calculate model min and max values over observation time
 * Simulates neuron model with finest dt and analyzes output range
 */
template <typename NeuronType, typename NeuronParamsType, typename StateType>
MinMaxResult calculate_min_max(
    NeuronModel model,
    const NeuronParamsType &params,
    const StateType &initial_state,
    double observation_time,
    double dt,
    double stabilization_time)
{
    // Validate time parameters
    if (observation_time <= 0 || dt <= 0 || stabilization_time < 0)
    {
        throw std::runtime_error("observation_time and dt must be positive, stabilization_time non-negative");
    }

    // Configure neuron model with provided parameters
    typename NeuronType::ConstructorArgs args;
    if (model == HINDMARSH_ROSE)
    {
        args.params[NeuronType::e] = params.e;
        args.params[NeuronType::mu] = params.mu;
        args.params[NeuronType::S] = params.S;
        args.params[NeuronType::a] = params.a;
        args.params[NeuronType::b] = params.b;
        args.params[NeuronType::c] = params.c;
        args.params[NeuronType::d] = params.d;
        args.params[NeuronType::xr] = params.xr;
        args.params[NeuronType::vh] = params.vh;
    }
    else
    {
        throw std::runtime_error("Model not implemented.");
    }

    NeuronType neuron(args);

    // Initialize neuron to specified state
    if (model == HINDMARSH_ROSE)
    {
        neuron.set(NeuronType::x, initial_state.x);
        neuron.set(NeuronType::y, initial_state.y);
        neuron.set(NeuronType::z, initial_state.z);
    }

    // Run stabilization phase to remove transients
    size_t stabilization_steps = static_cast<size_t>(stabilization_time / dt);
    for (size_t i = 0; i < stabilization_steps; i++)
        neuron.step(dt);

    // Find absolute min/max over observation time using finest dt
    double min = SignalConstants::DOUBLE_MAX;
    double max = SignalConstants::DOUBLE_MIN;
    size_t obs_steps = static_cast<size_t>(observation_time / dt);
    for (size_t step = 0; step < obs_steps; step++)
    {
        neuron.step(dt);

        double val;
        if (model == HINDMARSH_ROSE)
        {
            val = neuron.get(NeuronType::x);
        }

        if (val > max)
            max = val;
        if (val < min)
            min = val;
    }

    MinMaxResult result;
    result.min = min;
    result.max = max;

    return result;
}

/*
 * Calculate points per burst for each dt given min and max
 * Simulates neuron model with different time steps and analyzes bursting behavior
 */
template <typename NeuronType, typename NeuronParamsType, typename StateType, size_t N>
PtsResult calculate_pts(
    NeuronModel model,
    const NeuronParamsType &params,
    const StateType &initial_state,
    const std::array<double, N> &dts,
    double observation_time,
    double stabilization_time,
    double min_val,
    double max_val)
{
    // Configure neuron model with provided parameters
    typename NeuronType::ConstructorArgs args;
    if (model == HINDMARSH_ROSE)
    {
        args.params[NeuronType::e] = params.e;
        args.params[NeuronType::mu] = params.mu;
        args.params[NeuronType::S] = params.S;
        args.params[NeuronType::a] = params.a;
        args.params[NeuronType::b] = params.b;
        args.params[NeuronType::c] = params.c;
        args.params[NeuronType::d] = params.d;
        args.params[NeuronType::xr] = params.xr;
        args.params[NeuronType::vh] = params.vh;
    }
    else
    {
        throw std::runtime_error("Model not implemented.");
    }

    NeuronType neuron(args);

    PtsResult result;

    // Compute relative thresholds for burst detection (10% and 90% of range)
    double range = max_val - min_val;
    double th_on = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range + min_val;
    double th_up = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range + min_val;

    // Calculate points per burst for each dt in the array
    for (size_t i = 0; i < N; i++)
    {
        double dt = dts[i];

        // Reset neuron to initial state for this dt
        if (model == HINDMARSH_ROSE)
        {
            neuron.set(NeuronType::x, initial_state.x);
            neuron.set(NeuronType::y, initial_state.y);
            neuron.set(NeuronType::z, initial_state.z);
        }

        // Stabilization phase with current dt
        size_t stabilization_steps = static_cast<size_t>(stabilization_time / dts[i]);
        for (size_t j = 0; j < stabilization_steps; j++)
            neuron.step(dt);

        // Detect bursts using threshold crossings
        bool up;
        if (model == HINDMARSH_ROSE)
        {
            up = (neuron.get(NeuronType::x) > th_up);
        }
        double total_steps = 0;
        int bursts_seen = -1; // Start at -1 to skip partial first burst
        size_t steps_in_current_burst = 0;
        double act_time = 0.0;

        // Count steps per burst over observation period
        while (bursts_seen < ConstCalculatorConstants::BURSTS_TO_AVERAGE && act_time < observation_time)
        {
            neuron.step(dt);
            act_time += dt;

            double val;
            if (model == HINDMARSH_ROSE)
            {
                val = neuron.get(NeuronType::x);
            }

            // Detect burst onset (rising edge)
            if (!up && val > th_up)
            {
                up = true;
                bursts_seen++;
                total_steps += steps_in_current_burst;
                steps_in_current_burst = 0;
            }
            // Detect burst end (falling edge)
            else if (up && val < th_on)
            {
                up = false;
            }

            steps_in_current_burst++;
        }

        // Store average points per burst, or mark as invalid
        if (bursts_seen <= 0)
        {
            result.pts[i] = SignalConstants::DOUBLE_MAX; // Sentinel for invalid
            result.invalid_dts.push_back(dts[i]);
        }
        else
        {
            result.pts[i] = total_steps / static_cast<double>(bursts_seen);
        }
    }

    return result;
}

/*
 * Convert NumericIntegrator enum to string
 */
std::string integrator_to_string(NumericIntegrator integrator)
{
    switch (integrator)
    {
    case RK4:
        return "RK4";
    default:
        return "UNKNOWN";
    }
}

/*
 * Print DTS array, PTS array, and invalid dts from PtsResult
 */
void print_tables(const PtsResult &pr, NumericIntegrator integrator)
{
    std::string integrator_str = integrator_to_string(integrator);

    // Output DTS array with formatting
    size_t ds_size_minus_1 = DTS_SIZE - 1;
    std::cout << "inline constexpr std::array<double, " << DTS_SIZE << "> DTS_" << integrator_str << " = {";
    for (size_t i = 0; i < DTS_SIZE; i++)
    {
        if (i % 8 == 0)
            std::cout << "\n    ";
        std::cout << ConstCalculatorConfig::DTS[i];
        if (i < ds_size_minus_1)
            std::cout << ", ";
    }
    std::cout << "};\n";

    // Output PTS array with formatting
    std::cout << "inline constexpr std::array<double, " << DTS_SIZE << "> PTS_" << integrator_str << " = {";
    for (size_t i = 0; i < DTS_SIZE; i++)
    {
        if (i % 8 == 0)
            std::cout << "\n    ";
        std::cout << pr.pts[i];
        if (i < ds_size_minus_1)
            std::cout << ", ";
    }
    std::cout << "};\n";

    // Report any dts where burst detection failed
    if (!pr.invalid_dts.empty())
    {
        std::cout << "Invalid dts (no bursts detected): ";
        for (double dt : pr.invalid_dts)
        {
            std::cout << dt << " ";
        }
        std::cout << "\n";
    }
}

int main()
{
    // Calculate min and max for Hindmarsh-Rose model with RK4 integration
    MinMaxResult mmr = calculate_min_max<DifferentialNeuronWrapper<SystemWrapper<HindmarshRoseModel<double>>, RungeKutta4>, HindmarshRoseParams, HindmarshRoseState>(
        NeuronModel::HINDMARSH_ROSE,
        ConstCalculatorConfig::HR_PARAMS,
        ConstCalculatorConfig::HR_INITIAL_STATE,
        ConstCalculatorConfig::OBSERVATION_TIME,
        ConstCalculatorConfig::MINMAX_DT,
        ConstCalculatorConfig::STABILIZATION_TIME);

    // Output results in C++ format for direct inclusion in code
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "inline constexpr double MIN = " << mmr.min << ";\n";
    std::cout << "inline constexpr double MAX = " << mmr.max << ";\n";

    // Calculate pts for Hindmarsh-Rose model with RK4 integration
    PtsResult pr = calculate_pts<DifferentialNeuronWrapper<SystemWrapper<HindmarshRoseModel<double>>, RungeKutta4>, HindmarshRoseParams, HindmarshRoseState>(
        NeuronModel::HINDMARSH_ROSE,
        ConstCalculatorConfig::HR_PARAMS,
        ConstCalculatorConfig::HR_INITIAL_STATE,
        ConstCalculatorConfig::DTS,
        ConstCalculatorConfig::OBSERVATION_TIME,
        ConstCalculatorConfig::STABILIZATION_TIME,
        mmr.min,
        mmr.max);

    // Print DTS, PTS, and invalid dts using the new function
    print_tables(pr, NumericIntegrator::RK4);

    return 0;
}