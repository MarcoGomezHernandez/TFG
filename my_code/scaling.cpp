// Signal scaling and transformation for neural model integration
#include "scaling.h"
#include "scaling_utils.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <array>
#include <cstddef>
#include <stdexcept>

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
 * Hindmarsh-Rose model constants and precomputed lookup tables
 * MIN/MAX: Range of model output
 * DTS: Available time steps for integration
 * PTS: Points per burst for each corresponding dt
 */
struct HindmarshRose
{
    static constexpr double MIN = -1.668473;
    static constexpr double MAX = 1.764310;
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
    static constexpr std::array<double, 144> PTS = {
        1689645.000000, 1408038.000000, 1206890.000000, 1056028.000000, 938692.000000, 844822.500000, 768021.000000, 704019.000000,
        649863.500000, 603445.000000, 563215.000000, 528014.000000, 469346.000000, 422411.000000, 384010.500000, 337929.000000,
        301722.500000, 291318.000000, 281607.500000, 272523.500000, 264007.000000, 256007.000000, 248477.000000, 241378.000000,
        234673.000000, 228330.500000, 222322.000000, 216621.000000, 211205.500000, 206054.500000, 201148.500000, 196470.500000,
        192005.000000, 187738.500000, 183657.000000, 179749.500000, 176005.000000, 172412.500000, 168964.500000, 165651.500000,
        162466.000000, 156448.500000, 150861.500000, 145659.000000, 140803.500000, 136262.000000, 132003.500000, 128003.500000,
        124238.500000, 120689.000000, 117336.500000, 114165.000000, 109717.000000, 105602.500000, 101786.000000, 98235.000000,
        94923.500000, 91828.500000, 88002.500000, 84482.000000, 81233.000000, 77506.500000, 74107.500000, 70993.500000,
        67585.500000, 64490.500000, 61219.000000, 57864.500000, 54858.500000, 51829.500000, 48834.000000, 45666.000000,
        42453.500000, 39294.000000, 36258.500000, 33130.000000, 30064.500000, 29747.000000, 29436.000000, 29132.000000,
        28735.500000, 28349.500000, 27974.000000, 27608.500000, 27252.000000, 26905.000000, 26566.500000, 26237.000000,
        25915.000000, 25600.500000, 25294.000000, 24920.500000, 24558.500000, 24207.000000, 23865.000000, 23532.500000,
        23209.000000, 22895.000000, 22589.000000, 22232.000000, 21886.500000, 21551.500000, 21226.500000, 20911.500000,
        20605.000000, 20259.500000, 19925.000000, 19601.500000, 19287.500000, 18984.500000, 18649.000000, 18326.000000,
        18013.000000, 17711.000000, 17382.500000, 17066.500000, 16762.000000, 16436.000000, 16122.500000, 15820.000000,
        15500.500000, 15194.000000, 14873.000000, 14565.500000, 14246.000000, 13940.000000, 13625.000000, 13324.500000,
        13016.500000, 12703.000000, 12386.500000, 12085.000000, 11781.500000, 11477.500000, 11173.000000, 10871.000000,
        10571.500000, 10263.000000, 9960.500000, 9652.500000, 9353.000000, 9052.000000, 8751.500000, 8444.500000};
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
            std::string value;
            size_t current_col = 0;

            // Extract specified column
            while (std::getline(ss, value, ','))
            {
                if (current_col == column_index)
                {
                    data.push_back(std::stod(value));
                    break;
                }
                current_col++;
            }
        }

        current_line++;
        if (current_line >= end_index)
            break;
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
        double value = signal[i];
        if (!up && value > th_up)
        {
            changes++;
            up = true;
        }
        else if (up && value < th_on)
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
template <size_t N>
DTSelection set_pts_burst(const std::array<double, N> &dts,
                          const std::array<double, N> &pts,
                          double pts_live,
                          NumericIntegrator method)
{
    DTSelection selection;
    if (method == NumericIntegrator::RK4)
    {
        selection = select_dt_neuron_model(dts, pts, pts_live);
    }
    else
    {
        throw std::runtime_error("Unsupported integrator");
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
        selection = set_pts_burst(HindmarshRose::DTS, HindmarshRose::PTS, external_pts_per_burst, integrator);
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
        double max_window = SignalConstants::DOUBLE_MIN;
        double min_window = SignalConstants::DOUBLE_MAX;
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
                max_window != SignalConstants::DOUBLE_MIN && min_window != SignalConstants::DOUBLE_MAX)
            {
                drift_counter = 0;

                // Update scaling based on observed drift
                factors = fix_drift(min_abs_model, max_abs_model, min_window, max_window, stats);

                // Reset window trackers
                max_window = SignalConstants::DOUBLE_MIN;
                min_window = SignalConstants::DOUBLE_MAX;
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

    // Apply horizontal scaling (time interpolation)
    // Calculate output size: between each pair of points, insert (s_points-1) interpolated points
    size_t interpolated_size = ((signal_size - 1) * s_points) + 1;
    std::vector<double> interpolated_signal;
    interpolated_signal.reserve(interpolated_size);

    for (size_t i = 0; i < signal_size - 1; i++)
    {
        // Add original point
        interpolated_signal.push_back(signal[i]);

        // Add linearly interpolated points between current and next
        for (double j = 1.0; j < s_points; j++)
        {
            double alpha = j / s_points;
            double interp_val = signal[i] + (alpha * (signal[i + 1] - signal[i]));
            interpolated_signal.push_back(interp_val);
        }
    }

    // Add final point
    interpolated_signal.push_back(signal.back());

    result.scaled_signal = interpolated_signal;
    result.success = true;

    return result;
}
