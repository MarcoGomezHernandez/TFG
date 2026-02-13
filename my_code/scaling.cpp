// Signal scaling and transformation for neural model integration
#include "scaling.hpp"

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
 * Fills the provided vector with data points from specified column
 */
void read_csv_column(std::vector<double> &data, const std::string &csv_path, size_t column_index,
                     size_t start_index, size_t num_points)
{
    std::ifstream file(csv_path);

    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open CSV file: " + csv_path);
    }

    data.reserve(num_points);

    const size_t end_index = start_index + num_points;

    std::string line;
    std::string val;
    std::stringstream ss;
    size_t current_line = 0;

    // Parse CSV and extract target column
    while (current_line < end_index && std::getline(file, line))
    {
        // Skip lines before start_index
        if (current_line >= start_index)
        {
            ss.str(line);
            ss.clear();
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
    }

    // Precalculate data size for efficiency
    const size_t data_size = data.size();

    // Check if fewer points were read than expected and warn
    if (data_size < num_points)
    {
        data.shrink_to_fit();
        std::cout << "Warning: Fewer data points read (" << data_size << ") than expected (" << num_points << ") from CSV file: " << csv_path << std::endl;
    }

    file.close();
}

/*
 * Calculate signal period by counting threshold crossings
 * Uses hysteresis (th_up for rising edge, th_on for falling edge)
 * Returns period in same units as tiempo_observacion
 */
double signal_period(double tiempo_observacion, const std::vector<double> &signal, size_t obs_points,
                     double th_up, double th_on)
{
    // Initial state based on first sample
    bool up = (signal[0] > th_up);
    double changes = 0.0;

    // Count upward threshold crossings (burst onsets)
    for (size_t i = 0; i < obs_points; i++)
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
    const double rg_virtual = max_virtual - min_virtual;
    const double rg_viva = max_viva - min_viva;

    ScalingFactors factors;
    // Calculate slope
    const double scale_real_to_virtual = rg_virtual / rg_viva;
    factors.scale_real_to_virtual = scale_real_to_virtual; // Calculate offset to align minima
    // Calculate y-intercept
    factors.offset_real_to_virtual = min_virtual - (min_viva * scale_real_to_virtual);

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

    constexpr double INVALID_DT = SignalConstants::INVALID_DT;

    double dt_candidate = INVALID_DT;
    double pts_burst_candidate = SignalConstants::INVALID_PTS;

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
                dt_candidate = dts[i];
                pts_burst_candidate = pts[i];

                // Check if pts_burst is close to integer multiple of pts_live
                fractpart = std::modf(pts_burst_candidate / pts_live, &intpart);

                if (fractpart <= SignalPrivateConfig::DT_SELECTION_TOLERANCE * intpart)
                {
                    flag = true;
                }
                break;
            }
        }

        if (flag)
        {
            break;
        }
    }

    // If no good match found, use last candidate
    if (!flag)
    {
        for (size_t i = N - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {
                dt_candidate = dts[i];
                pts_burst_candidate = pts[i];
                break;
            }
        }
    }

    DTSelection selection;

    selection.success = (dt_candidate != INVALID_DT);
    selection.dt = dt_candidate;
    selection.pts_burst = pts_burst_candidate;

    return selection;
}

/*
 * Dispatch dt selection based on integration method
 * Currently only supports RK4
 */
DTSelection set_pts_burst(NeuronModel model, NumericIntegrator integrator, double pts_live)
{
    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        if (integrator == NumericIntegrator::RK4)
        {
            return select_dt_neuron_model(HindmarshRose::DTS_RK4, HindmarshRose::PTS_RK4, pts_live);
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

    constexpr double DRIFT_PERCENTAGE_MIN = SignalPrivateConfig::DRIFT_PERCENTAGE_MIN;
    constexpr double DRIFT_PERCENTAGE_MAX = SignalPrivateConfig::DRIFT_PERCENTAGE_MAX;

    // Adjust relative thresholds with drift margins
    // Add margin for positive values, subtract for negative
    if (min_window > 0)
    {
        stats.min_rel_real = min_window + (min_window * DRIFT_PERCENTAGE_MIN);
    }
    else
    {
        stats.min_rel_real = min_window - (min_window * DRIFT_PERCENTAGE_MIN);
    }

    if (max_window > 0)
    {
        stats.max_rel_real = max_window - (max_window * DRIFT_PERCENTAGE_MAX);
    }
    else
    {
        stats.max_rel_real = max_window + (max_window * DRIFT_PERCENTAGE_MAX);
    }

    return factors;
}

/*
 * Initialize signal statistics from observation window
 * Computes min/max, relative thresholds, and signal period
 */
SignalStats ini_recibido(const std::vector<double> &signal, size_t obs_points, double csv_step)
{
    const double observation_time_to_use = obs_points * csv_step;

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
    const double range = max_abs - min_abs;
    const double min_rel_real = SignalPublicConfig::SIGNAL_PERCENTAGE_MIN * range + min_abs;
    const double max_rel_real = SignalPublicConfig::SIGNAL_PERCENTAGE_MAX * range + min_abs;
    stats.min_rel_real = min_rel_real;
    stats.max_rel_real = max_rel_real;

    // Calculate signal period using threshold crossings
    stats.period_signal = signal_period(observation_time_to_use, signal, obs_points, max_rel_real, min_rel_real);

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
    // Validate all input parameters
    if (csv_step <= 0 || use_time <= 0 || observation_time <= 0 || start_time < 0 || column_index < 0 || csv_path.empty())
    {
        throw std::runtime_error("Invalid arguments: csv_step, use_time, observation_time must be positive, start_time and column_index non-negative, csv_path non-empty");
    }

    ScaledSignalResult result;

    // References to result members to avoid direct struct access
    std::vector<double> &signal = result.signal;
    std::vector<double> &interpolated_points = result.interpolated_points;

    // Calculate indices and points once
    const size_t start_index = static_cast<size_t>(start_time / csv_step);
    const size_t use_points = static_cast<size_t>(use_time / csv_step);
    if (use_points == 0)
    {
        throw std::runtime_error("use_time is too short to read any points with given csv_step");
    }

    const size_t obs_points = static_cast<size_t>(observation_time / csv_step);
    if (obs_points == 0)
    {
        throw std::runtime_error("observation_time is too short to read any points with given csv_step");
    }
    const size_t read_points = std::max(use_points, obs_points);

    // Read signal data from CSV file directly into result.signal
    read_csv_column(signal, csv_path, column_index, start_index, read_points);

    size_t signal_size = signal.size();
    if (signal_size == 0)
    {
        throw std::runtime_error("No data read from CSV file");
    }

    // Extract signal statistics from observation window
    SignalStats stats = ini_recibido(signal, obs_points, csv_step);

    // Trim signal to the size for use_time after observation
    if (signal_size > use_points)
    {
        signal.resize(use_points);
        signal_size = use_points;
    }

    // Calculate points per burst in external signal
    const double external_pts_per_burst = stats.period_signal / csv_step;

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
        signal.clear();
        result.dt = SignalConstants::INVALID_DT;
        return result;
    }

    result.dt = selection.dt;

    // Calculate horizontal scaling factor (time interpolation)
    size_t s_points = static_cast<size_t>(selection.pts_burst / external_pts_per_burst);
    if (s_points == 0)
        s_points = 1;

    result.points_factor = s_points;

    double &min_abs_real = stats.min_abs_real;
    double &max_abs_real = stats.max_abs_real;

    // Calculate initial vertical scaling factors
    ScalingFactors factors = calcula_escala(min_abs_model, max_abs_model, min_abs_real, max_abs_real);
    double scale_real_to_virtual = factors.scale_real_to_virtual;
    double offset_real_to_virtual = factors.offset_real_to_virtual;

    // Apply vertical scaling (amplitude transformation)
    if (check_drift)
    {
        // Drift checking mode: recalculate scaling periodically
        size_t drift_counter = 0;
        constexpr double DOUBLE_MAX = GeneralConstants::DOUBLE_MAX;
        constexpr double DOUBLE_MIN = GeneralConstants::DOUBLE_MIN;
        double max_window = DOUBLE_MIN;
        double min_window = DOUBLE_MAX;
        const double drift_aux_range = max_abs_real - min_abs_real;

        for (size_t i = 0; i < signal_size; i++)
        {
            double val = signal[i];

            // Track windowed min/max within reasonable bounds
            if ((min_window > val) && (val > (min_abs_real - drift_aux_range)))
            {
                min_window = val;
            }
            if ((max_window < val) && (val < (max_abs_real + drift_aux_range)))
            {
                max_window = val;
            }

            // Recalculate scaling factors every N bursts
            if (drift_counter >= (SignalPrivateConfig::DRIFT_N_BURST * external_pts_per_burst) &&
                max_window != DOUBLE_MIN && min_window != DOUBLE_MAX)
            {
                drift_counter = 0;

                // Update scaling based on observed drift
                factors = fix_drift(min_abs_model, max_abs_model, min_window, max_window, stats);
                scale_real_to_virtual = factors.scale_real_to_virtual;
                offset_real_to_virtual = factors.offset_real_to_virtual;

                // Reset window trackers
                max_window = DOUBLE_MIN;
                min_window = DOUBLE_MAX;
            }

            drift_counter++;

            // Apply vertical scaling transformation
            signal[i] = val * scale_real_to_virtual + offset_real_to_virtual;
        }
    }
    else
    {
        // Simple mode: apply constant scaling to all points
        for (size_t i = 0; i < signal_size; i++)
        {
            signal[i] = signal[i] * scale_real_to_virtual + offset_real_to_virtual;
        }
    }

    // Apply horizontal scaling (time interpolation)
    // Calculate output size: only the newly interpolated points between originals
    const size_t interpolated_size = (signal_size - 1) * (s_points - 1);
    interpolated_points.reserve(interpolated_size);

    for (size_t i = 0; i < signal_size - 1; i++)
    {
        // Add linearly interpolated points between current and next, excluding originals
        for (double j = 1.0; j < s_points; j++)
        {
            double alpha = j / s_points;
            double interp_val = signal[i] + (alpha * (signal[i + 1] - signal[i]));
            interpolated_points.push_back(interp_val);
        }
    }

    result.success = true;

    return result;
}
