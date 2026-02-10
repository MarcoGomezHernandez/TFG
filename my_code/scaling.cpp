// Signal scaling and transformation for neural model integration
#include "scaling.h"
#include "utils.h"
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

/*
 * Private configuration constants for signal processing algorithms
 */
namespace SignalPrivateConfig
{
    // Tolerance for dt selection algorithm
    static constexpr double DT_SELECTION_TOLERANCE = 0.1;

    // Drift correction margins (percentage of range)
    static constexpr double DRIFT_PERCENTAGE_MIN = 0.1;
    static constexpr double DRIFT_PERCENTAGE_MAX = 0.1;
    // Number of bursts before recalculating scaling factors
    static constexpr size_t DRIFT_N_BURST = 2;
}

/*
 * Scaling factors for linear transformation: y = scale * x + offset
 */
struct ScalingFactors
{
    double scale_real_to_virtual;  // Multiplicative scale factor
    double offset_real_to_virtual; // Additive offset
};

/*
 * Result of dt selection from lookup table
 */
struct DTSelection
{
    double dt;        // Selected time step
    double pts_burst; // Points per burst for selected dt
    bool success;     // Whether a valid dt was found
};

/*
 * Statistical properties extracted from input signal
 */
struct SignalStats
{
    double min_abs_real;  // Absolute minimum value
    double max_abs_real;  // Absolute maximum value
    double min_rel_real;  // Relative minimum threshold (10% of range)
    double max_rel_real;  // Relative maximum threshold (90% of range)
    double period_signal; // Calculated signal period
};

/*
 * Read specific column from CSV file within time range
 * Returns vector containing data points from specified column
 */
std::vector<double> read_csv_column(const std::string &csv_path, size_t column_index,
                                    double start_time, double use_time, double csv_step)
{
    std::vector<double> data;
    std::ifstream file(csv_path);

    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open CSV file: " + csv_path);
    }

    // Convert time parameters to array indices
    size_t start_index = static_cast<size_t>(start_time / csv_step);
    size_t num_points = static_cast<size_t>(use_time / csv_step);

    data.reserve(num_points);

    size_t end_index = start_index + num_points;

    std::string line;
    size_t current_line = 0;

    // Parse CSV and extract target column
    while (std::getline(file, line))
    {
        // Skip lines before start_index
        if (current_line >= start_index)
        {
            std::stringstream ss(line);
            std::string val;
            size_t current_col = 0;

            // Extract specified column
            while (std::getline(ss, val, ','))
            {
                if (current_col == column_index)
                {
                    data.push_back(std::stod(val));
                    break;
                }
                current_col++;
            }
        }

        current_line++;
        if (current_line >= end_index)
            break;
    }

    // Precalculate data size for efficiency
    size_t data_size = data.size();

    // Check if fewer points were read than expected and warn
    if (data_size < num_points)
    {
        std::cout << "Warning: Fewer data points read (" << data_size << ") than expected (" << num_points << ") from CSV file: " << csv_path << std::endl;
    }

    file.close();
    return data;
}

/*
 * Calculate signal period by counting threshold crossings
 * Uses hysteresis (th_up for rising edge, th_on for falling edge)
 * Returns period in same units as tiempo_observacion
 */
double signal_period(double tiempo_observacion, const std::vector<double> &signal, size_t size,
                     double th_up, double th_on)
{
    // Initial state based on first sample
    bool up = (signal[0] > th_up);
    double changes = 0.0;

    // Count upward threshold crossings (burst onsets)
    for (size_t i = 0; i < size; i++)
    {
        double val = signal[i];
        if (!up && val > th_up)
        {
            changes++;
            up = true;
        }
        else if (up && val < th_on)
        {
            up = false;
        }
    }

    // Period = observation_time / number_of_bursts
    return 1.0 / (changes / tiempo_observacion);
}

/*
 * Calculate linear scaling factors to map one range to another
 * Maps [min_viva, max_viva] -> [min_virtual, max_virtual]
 */
ScalingFactors calcula_escala(double min_virtual, double max_virtual,
                              double min_viva, double max_viva)
{
    double rg_virtual = max_virtual - min_virtual;
    double rg_viva = max_viva - min_viva;

    ScalingFactors factors;
    // Calculate slope
    factors.scale_real_to_virtual = rg_virtual / rg_viva;
    // Calculate y-intercept
    factors.offset_real_to_virtual = min_virtual - (min_viva * factors.scale_real_to_virtual);

    return factors;
}

/*
 * Select appropriate dt from lookup table matching live signal period
 * Searches for dt where pts_burst is close to a multiple of pts_live
 * Uses tolerance to find acceptable matches
 */
template <size_t N>
DTSelection select_dt_neuron_model(const std::array<double, N> &dts,
                                   const std::array<double, N> &pts,
                                   double pts_live)
{
    double aux = pts_live;
    double factor = 1.0;
    double intpart, fractpart;
    bool flag = false;

    DTSelection selection;
    selection.dt = SignalConstants::INVALID_DT;
    selection.pts_burst = SignalConstants::INVALID_PTS;

    // Search for matching dt by scaling pts_live
    while (aux < pts[0])
    {
        aux = pts_live * factor;
        factor += 1.0;

        // Search backwards through lookup table
        for (size_t i = N - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {
                selection.dt = dts[i];
                selection.pts_burst = pts[i];

                // Check if pts_burst is close to integer multiple of pts_live
                fractpart = std::modf(selection.pts_burst / pts_live, &intpart);

                if (fractpart <= SignalPrivateConfig::DT_SELECTION_TOLERANCE * intpart)
                {
                    flag = true;
                }
                break;
            }
        }

        if (flag)
            break;
    }

    // If no good match found, use last candidate
    if (!flag)
    {
        for (size_t i = N - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {
                selection.dt = dts[i];
                selection.pts_burst = pts[i];
                break;
            }
        }
    }

    selection.success = (selection.dt != SignalConstants::INVALID_DT);

    return selection;
}

/*
 * Dispatch dt selection based on integration method
 * Currently only supports RK4
 */
DTSelection set_pts_burst(NeuronModel model, NumericIntegrator integrator, double pts_live)
{
    DTSelection selection;
    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        if (integrator == NumericIntegrator::RK4)
        {
            selection = select_dt_neuron_model(HindmarshRose::DTS_RK4, HindmarshRose::PTS_RK4, pts_live);
        }
        else
        {
            throw std::runtime_error("Unsupported integrator");
        }
    }
    else
    {
        throw std::runtime_error("Unsupported neuron model");
    }
    return selection;
}

/*
 * Adjust scaling factors to prevent signal drift
 * Recalculates scaling based on recent min/max window
 * Updates relative thresholds with safety margins
 */
ScalingFactors fix_drift(double min_abs_model, double max_abs_model, double min_window, double max_window, SignalStats &stats)
{
    // Recalculate scaling based on observed window
    ScalingFactors factors = calcula_escala(min_abs_model, max_abs_model, min_window, max_window);

    // Adjust relative thresholds with drift margins
    // Add margin for positive values, subtract for negative
    if (min_window > 0)
    {
        stats.min_rel_real = min_window + (min_window * SignalPrivateConfig::DRIFT_PERCENTAGE_MIN);
    }
    else
    {
        stats.min_rel_real = min_window - (min_window * SignalPrivateConfig::DRIFT_PERCENTAGE_MIN);
    }

    if (max_window > 0)
    {
        stats.max_rel_real = max_window - (max_window * SignalPrivateConfig::DRIFT_PERCENTAGE_MAX);
    }
    else
    {
        stats.max_rel_real = max_window + (max_window * SignalPrivateConfig::DRIFT_PERCENTAGE_MAX);
    }

    return factors;
}

/*
 * Initialize signal statistics from observation window
 * Computes min/max, relative thresholds, and signal period
 */
SignalStats ini_recibido(const std::vector<double> &signal, double observation_time, double csv_step)
{
    size_t signal_size = signal.size();

    // Limit observation to available data
    size_t obs_points = static_cast<size_t>(observation_time / csv_step);
    if (obs_points > signal_size)
        obs_points = signal_size;

    double observation_time_to_use = obs_points * csv_step;

    SignalStats stats;

    // Find absolute min/max in observation window
    double max_abs = -DBL_MAX;
    double min_abs = DBL_MAX;

    for (size_t i = 0; i < obs_points; i++)
    {
        double val = signal[i];
        if (val > max_abs)
            max_abs = val;
        if (val < min_abs)
            min_abs = val;
    }

    stats.min_abs_real = min_abs;
    stats.max_abs_real = max_abs;

    // Calculate relative thresholds (10% and 90% of range)
    double range = max_abs - min_abs;
    stats.min_rel_real = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range + min_abs;
    stats.max_rel_real = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range + min_abs;

    // Calculate signal period using threshold crossings
    stats.period_signal = signal_period(observation_time_to_use, signal, obs_points, stats.max_rel_real, stats.min_rel_real);

    return stats;
}

/*
 * Main function: scale external signal to match neural model dynamics
 * Performs both vertical scaling (amplitude) and horizontal scaling (time)
 * Optional drift correction recalculates scaling factors periodically
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
    bool check_drift)
{
    ScaledSignalResult result;

    // Validate all input parameters
    if (csv_step <= 0 || use_time <= 0 || observation_time <= 0 || start_time < 0 || column_index < 0 || csv_path.empty())
    {
        throw std::runtime_error("Invalid arguments: csv_step, use_time, observation_time must be positive, start_time and column_index non-negative, csv_path non-empty");
    }

    // Read signal data from CSV file
    std::vector<double> signal = read_csv_column(csv_path, column_index, start_time, use_time, csv_step);

    size_t signal_size = signal.size();
    if (signal_size == 0)
    {
        throw std::runtime_error("No data read from CSV file");
    }

    // Extract signal statistics from observation window
    SignalStats stats = ini_recibido(signal, observation_time, csv_step);

    // Calculate points per burst in external signal
    double external_pts_per_burst = stats.period_signal / csv_step;

    // Select appropriate dt and get model parameters
    DTSelection selection;
    double min_abs_model, max_abs_model;
    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        min_abs_model = HindmarshRose::MIN;
        max_abs_model = HindmarshRose::MAX;
        selection = set_pts_burst(model, integrator, external_pts_per_burst);
    }
    else
    {
        throw std::runtime_error("Unsupported neuron model");
    }

    // Check if valid dt was found
    if (!selection.success)
    {
        result.success = false;
        result.dt = SignalConstants::INVALID_DT;
        return result;
    }

    result.dt = selection.dt;

    // Calculate horizontal scaling factor (time interpolation)
    size_t s_points = static_cast<size_t>(selection.pts_burst / external_pts_per_burst);
    if (s_points == 0)
        s_points = 1;

    // Calculate initial vertical scaling factors
    ScalingFactors factors = calcula_escala(min_abs_model, max_abs_model, stats.min_abs_real, stats.max_abs_real);

    // Apply vertical scaling (amplitude transformation)
    if (check_drift)
    {
        // Drift checking mode: recalculate scaling periodically
        size_t drift_counter = 0;
        double max_window = GeneralConstants::DOUBLE_MIN;
        double min_window = GeneralConstants::DOUBLE_MAX;
        double drift_aux_range = stats.max_abs_real - stats.min_abs_real;

        for (size_t i = 0; i < signal_size; i++)
        {
            double val = signal[i];

            // Track windowed min/max within reasonable bounds
            if ((min_window > val) && (val > (stats.min_abs_real - drift_aux_range)))
            {
                min_window = val;
            }
            if ((max_window < val) && (val < (stats.max_abs_real + drift_aux_range)))
            {
                max_window = val;
            }

            // Recalculate scaling factors every N bursts
            if (drift_counter >= (SignalPrivateConfig::DRIFT_N_BURST * external_pts_per_burst) &&
                max_window != GeneralConstants::DOUBLE_MIN && min_window != GeneralConstants::DOUBLE_MAX)
            {
                drift_counter = 0;

                // Update scaling based on observed drift
                factors = fix_drift(min_abs_model, max_abs_model, min_window, max_window, stats);

                // Reset window trackers
                max_window = GeneralConstants::DOUBLE_MIN;
                min_window = GeneralConstants::DOUBLE_MAX;
            }

            drift_counter++;

            // Apply vertical scaling transformation
            signal[i] = val * factors.scale_real_to_virtual + factors.offset_real_to_virtual;
        }
    }
    else
    {
        // Simple mode: apply constant scaling to all points
        for (size_t i = 0; i < signal_size; i++)
        {
            signal[i] = signal[i] * factors.scale_real_to_virtual + factors.offset_real_to_virtual;
        }
    }

    // Store non-interpolated scaled signal
    result.signal = signal;
    result.points_factor = s_points;

    // Apply horizontal scaling (time interpolation)
    // Calculate output size: only the newly interpolated points between originals
    size_t interpolated_size = (signal_size - 1) * (s_points - 1);
    std::vector<double> interpolated_signal;
    interpolated_signal.reserve(interpolated_size);

    for (size_t i = 0; i < signal_size - 1; i++)
    {
        // Add linearly interpolated points between current and next, excluding originals
        for (double j = 1.0; j < s_points; j++)
        {
            double alpha = j / s_points;
            double interp_val = signal[i] + (alpha * (signal[i + 1] - signal[i]));
            interpolated_signal.push_back(interp_val);
        }
    }

    result.interpolated_points = interpolated_signal;
    result.success = true;

    return result;
}
