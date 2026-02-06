#include "scaling.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <array>
#include <cstddef>

// Structs for return values
struct ScalingFactors
{
    double scale_real_to_virtual;
    double offset_real_to_virtual;
};

struct DTSelection
{
    double dt;
    double pts_burst;
    bool success;
};

struct SignalStats
{
    double min_abs_real;
    double max_abs_real;
    double min_rel_real;
    double max_rel_real;
    double period_signal;
};

struct HindmarshRose
{
    static constexpr double min = -1.608734;
    static constexpr double max = 1.797032;
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
    static constexpr std::array<double, 144> pts = {1861445.200000, 1551204.200000, 1329603.600000, 1163403.200000, 1034136.200000, 930722.600000, 846111.400000, 775602.200000, 715940.400000, 664801.800000, 620481.800000, 581701.800000, 517068.200000, 465361.200000, 423055.600000, 372289.000000, 332401.000000, 320938.800000, 310240.800000, 300233.000000, 290850.800000, 282037.200000, 273742.000000, 265920.800000, 258534.000000, 251546.600000, 244927.000000, 238646.800000, 232680.600000, 227005.400000, 221600.600000, 216447.000000, 211527.800000, 206827.400000, 202331.000000, 198026.200000, 193900.600000, 189943.200000, 186144.400000, 182494.600000, 178985.200000, 172356.000000, 166200.400000, 160469.400000, 155120.400000, 150116.600000, 145425.400000, 141018.400000, 136871.000000, 132960.400000, 129267.000000, 125773.200000, 120873.000000, 116340.400000, 112135.400000, 108223.600000, 104575.600000, 101165.600000, 96950.200000, 93072.200000, 89492.600000, 85387.400000, 81642.400000, 78212.000000, 74457.800000, 71047.600000, 67443.600000, 63748.200000, 60436.600000, 57099.600000, 53489.800000, 50038.800000, 46536.200000, 43089.000000, 39605.000000, 36215.000000, 32887.600000, 32542.800000, 32205.000000, 31765.200000, 31337.400000, 30921.000000, 30515.400000, 30120.400000, 29735.600000, 29360.400000, 28994.400000, 28637.400000, 28289.400000, 27949.400000, 27617.800000, 27214.000000, 26821.800000, 26440.800000, 26070.600000, 25710.600000, 25360.200000, 25019.200000, 24687.400000, 24300.800000, 23925.800000, 23562.600000, 23210.000000, 22867.800000, 22535.600000, 22160.000000, 21796.600000, 21445.200000, 21104.800000, 20728.800000, 20365.800000, 20015.400000, 19676.800000, 19309.400000, 18955.600000, 18614.200000, 18249.200000, 17898.400000, 17560.600000, 17203.400000, 16860.600000, 16502.000000, 16158.000000, 15801.400000, 15460.200000, 15108.800000, 14773.000000, 14429.400000, 14080.000000, 13747.400000, 13410.400000, 13071.400000, 12731.400000, 12392.400000, 12055.200000, 11706.200000, 11363.200000, 11026.400000, 10684.400000, 10340.200000, 9995.600000, 9653.200000, 9305.400000};
};

std::vector<double> read_csv_column(const std::string &csv_path, size_t column_index,
                                    double start_time, double use_time, double csv_step)
{
    std::vector<double> data;
    std::ifstream file(csv_path);

    if (!file.is_open())
    {
        return data;
    }

    size_t start_index = static_cast<size_t>(start_time / csv_step);
    size_t num_points = static_cast<size_t>(use_time / csv_step);

    data.reserve(num_points);

    size_t end_index = start_index + num_points;

    std::string line;
    size_t current_line = 0;

    while (std::getline(file, line))
    {
        if (current_line >= start_index)
        {
            std::stringstream ss(line);
            std::string value;
            size_t current_col = 0;

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

double signal_period(double tiempo_observacion, const std::vector<double> &signal, size_t size,
                     double th_up, double th_on)
{
    bool up = (signal[0] > th_up);
    double changes = 0.0;

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

    return 1.0 / (changes / tiempo_observacion);
}

ScalingFactors calcula_escala(double min_virtual, double max_virtual,
                              double min_viva, double max_viva)
{
    double rg_virtual = max_virtual - min_virtual;
    double rg_viva = max_viva - min_viva;

    ScalingFactors factors;
    factors.scale_real_to_virtual = rg_virtual / rg_viva;
    factors.offset_real_to_virtual = min_virtual - (min_viva * factors.scale_real_to_virtual);

    return factors;
}

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
    selection.dt = -1.0;
    selection.pts_burst = -1.0;

    while (aux < pts[0])
    {
        aux = pts_live * factor;
        factor += 1.0;

        for (size_t i = N - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {
                selection.dt = dts[i];
                selection.pts_burst = pts[i];

                fractpart = std::modf(selection.pts_burst / pts_live, &intpart);

                if (fractpart <= 0.1 * intpart)
                {
                    flag = true;
                }
                break;
            }
        }

        if (flag)
            break;
    }

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

    selection.success = (selection.dt != -1.0);

    return selection;
}

template <size_t N>
DTSelection set_pts_burst(const std::array<double, N> &dts,
                          const std::array<double, N> &pts,
                          double pts_live,
                          Integrator method)
{
    DTSelection selection;
    if (method == Integrator::RK4)
    {
        selection = select_dt_neuron_model(dts, pts, pts_live);
    }
    else
    {
        // Invalid for unsupported integrators
        selection.dt = -1.0;
        selection.pts_burst = -1.0;
        selection.success = false;
    }
    return selection;
}

ScalingFactors fix_drift(double min_abs_model, double max_abs_model, double min_window, double max_window, SignalStats &stats)
{
    ScalingFactors factors = calcula_escala(min_abs_model, max_abs_model, min_window, max_window);

    double per_min = 0.1, per_max = 0.1;
    if (min_window > 0)
    {
        stats.min_rel_real = min_window + (min_window * per_min);
    }
    else
    {
        stats.min_rel_real = min_window - (min_window * per_min);
    }

    if (max_window > 0)
    {
        stats.max_rel_real = max_window - (max_window * per_max);
    }
    else
    {
        stats.max_rel_real = max_window + (max_window * per_max);
    }

    return factors;
}

ScaledSignalResult scale_signal(
    const std::string &csv_path,
    size_t column_index,
    double csv_step,
    double start_time,
    double use_time,
    double observation_time,
    Integrator integrator,
    NeuronModel model,
    bool check_drift)
{
    ScaledSignalResult result;
    result.success = false;
    result.dt = -1.0;

    if (csv_step <= 0 || use_time <= 0 || observation_time <= 0 || start_time < 0 || column_index < 0 || csv_path.empty())
    {
        return result;
    }

    // Read CSV data
    std::vector<double> signal = read_csv_column(csv_path, column_index, start_time, use_time, csv_step);

    size_t signal_size = signal.size();
    if (signal_size == 0)
    {
        return result;
    }

    // Get signal statistics using ini_recibido
    SignalStats stats = ini_recibido(signal, observation_time, csv_step);

    // Calculate dt and pts_burst
    double external_pts_per_burst = stats.period_signal / csv_step;

    DTSelection selection;
    double min_abs_model, max_abs_model;
    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        min_abs_model = HindmarshRose::min;
        max_abs_model = HindmarshRose::max;
        selection = set_pts_burst(HindmarshRose::dts, HindmarshRose::pts, external_pts_per_burst, integrator);
    }
    else
    {
        return result;
    }

    if (!selection.success)
    {
        return result;
    }

    result.dt = selection.dt;

    // Calculate s_points
    size_t s_points = static_cast<size_t>(selection.pts_burst / external_pts_per_burst);
    if (s_points == 0)
        s_points = 1;

    // Initial scaling factors
    ScalingFactors factors = calcula_escala(min_abs_model, max_abs_model, stats.min_abs_real, stats.max_abs_real);

    // Apply vertical scaling in-place to signal
    if (check_drift)
    {
        // Drift checking logic with in-place scaling
        const size_t drift_n_burst = 2;
        size_t drift_counter = 0;
        double max_window = -999999.0;
        double min_window = 999999.0;
        double drift_aux_range = stats.max_abs_real - stats.min_abs_real;

        for (size_t i = 0; i < signal_size; i++)
        {
            double val = signal[i];
            // Update window
            if ((min_window > val) && (val > (stats.min_abs_real - drift_aux_range)))
            {
                min_window = val;
            }
            if ((max_window < val) && (val < (stats.max_abs_real + drift_aux_range)))
            {
                max_window = val;
            }

            // Recalculate every drift_n_burst bursts
            if (drift_counter >= (drift_n_burst * external_pts_per_burst) &&
                max_window != -999999.0 && min_window != 999999.0)
            {
                drift_counter = 0;

                factors = fix_drift(min_abs_model, max_abs_model, min_window, max_window, stats);

                max_window = -999999.0;
                min_window = 999999.0;
            }

            drift_counter++;

            // Scale current point in-place (vertical scaling)
            signal[i] = val * factors.scale_real_to_virtual + factors.offset_real_to_virtual;
        }
    }
    else
    {
        // No drift checking - simple scaling
        for (size_t i = 0; i < signal_size; i++)
        {
            signal[i] = signal[i] * factors.scale_real_to_virtual + factors.offset_real_to_virtual;
        }
    }

    // Perform linear interpolation (horizontal scaling) on the now-scaled signal
    size_t interpolated_size = ((signal_size - 1) * s_points) + 1;
    std::vector<double> interpolated_signal;
    interpolated_signal.reserve(interpolated_size);

    for (size_t i = 0; i < signal_size - 1; i++)
    {
        interpolated_signal.push_back(signal[i]);

        // Add s_points - 1 intermediate points
        for (double j = 1.0; j < s_points; j++)
        {
            double alpha = j / s_points;
            double interp_val = signal[i] + (alpha * (signal[i + 1] - signal[i]));
            interpolated_signal.push_back(interp_val);
        }
    }

    // Add last point
    interpolated_signal.push_back(signal.back());

    result.scaled_signal = interpolated_signal;
    result.success = true;

    return result;
}

SignalStats ini_recibido(const std::vector<double> &signal, double observation_time, double csv_step)
{
    size_t signal_size = signal.size();

    // Calculate observation range
    size_t obs_points = static_cast<size_t>(observation_time / csv_step);
    if (obs_points > signal_size)
        obs_points = signal_size;

    double observation_time_to_use = obs_points * csv_step;

    SignalStats stats;

    // Determine min/max from signal
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

    double range = max_abs - min_abs;
    double percentage_min = 0.10;
    double percentage_max = 0.90;
    stats.min_rel_real = percentage_min * range + min_abs;
    stats.max_rel_real = percentage_max * range + min_abs;

    // Calculate signal period
    stats.period_signal = signal_period(observation_time_to_use, signal, obs_points, stats.max_rel_real, stats.min_rel_real);

    return stats;
}
