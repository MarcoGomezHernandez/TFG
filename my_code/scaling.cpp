#include "scaling.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cfloat>

// Precomputed tables for Hindmarsh-Rose RK4
static const std::vector<double> HR_DTS_RK4 = {
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

static const std::vector<double> HR_PTS_RK4 = {1861445.200000, 1551204.200000, 1329603.600000, 1163403.200000, 1034136.200000, 930722.600000, 846111.400000, 775602.200000, 715940.400000, 664801.800000, 620481.800000, 581701.800000, 517068.200000, 465361.200000, 423055.600000, 372289.000000, 332401.000000, 320938.800000, 310240.800000, 300233.000000, 290850.800000, 282037.200000, 273742.000000, 265920.800000, 258534.000000, 251546.600000, 244927.000000, 238646.800000, 232680.600000, 227005.400000, 221600.600000, 216447.000000, 211527.800000, 206827.400000, 202331.000000, 198026.200000, 193900.600000, 189943.200000, 186144.400000, 182494.600000, 178985.200000, 172356.000000, 166200.400000, 160469.400000, 155120.400000, 150116.600000, 145425.400000, 141018.400000, 136871.000000, 132960.400000, 129267.000000, 125773.200000, 120873.000000, 116340.400000, 112135.400000, 108223.600000, 104575.600000, 101165.600000, 96950.200000, 93072.200000, 89492.600000, 85387.400000, 81642.400000, 78212.000000, 74457.800000, 71047.600000, 67443.600000, 63748.200000, 60436.600000, 57099.600000, 53489.800000, 50038.800000, 46536.200000, 43089.000000, 39605.000000, 36215.000000, 32887.600000, 32542.800000, 32205.000000, 31765.200000, 31337.400000, 30921.000000, 30515.400000, 30120.400000, 29735.600000, 29360.400000, 28994.400000, 28637.400000, 28289.400000, 27949.400000, 27617.800000, 27214.000000, 26821.800000, 26440.800000, 26070.600000, 25710.600000, 25360.200000, 25019.200000, 24687.400000, 24300.800000, 23925.800000, 23562.600000, 23210.000000, 22867.800000, 22535.600000, 22160.000000, 21796.600000, 21445.200000, 21104.800000, 20728.800000, 20365.800000, 20015.400000, 19676.800000, 19309.400000, 18955.600000, 18614.200000, 18249.200000, 17898.400000, 17560.600000, 17203.400000, 16860.600000, 16502.000000, 16158.000000, 15801.400000, 15460.200000, 15108.800000, 14773.000000, 14429.400000, 14080.000000, 13747.400000, 13410.400000, 13071.400000, 12731.400000, 12392.400000, 12055.200000, 11706.200000, 11363.200000, 11026.400000, 10684.400000, 10340.200000, 9995.600000, 9653.200000, 9305.400000};

std::vector<double> read_csv_column(const std::string &csv_path, int column_index,
                                    double start_time, double use_time, double csv_step)
{
    std::vector<double> data;
    std::ifstream file(csv_path);

    if (!file.is_open())
    {
        return data;
    }

    int start_index = static_cast<int>(start_time / csv_step);
    int num_points = static_cast<int>(use_time / csv_step);

    std::string line;
    int current_line = 0;

    while (std::getline(file, line))
    {
        if (current_line >= start_index && current_line < start_index + num_points)
        {
            std::stringstream ss(line);
            std::string value;
            int col = 0;

            while (std::getline(ss, value, ','))
            {
                if (col == column_index)
                {
                    data.push_back(std::stod(value));
                    break;
                }
                col++;
            }
        }

        current_line++;
        if (current_line >= start_index + num_points)
            break;
    }

    file.close();
    return data;
}

double signal_period(int seg_observacion, const std::vector<double> &signal,
                     double th_up, double th_on)
{
    bool up = (signal[0] > th_up);
    double changes = 0.0;

    for (const auto& value : signal)
    {
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

    return 1.0 / (changes / seg_observacion);
}

void calcula_escala(double min_virtual, double max_virtual,
                    double min_viva, double max_viva,
                    double &scale_virtual_to_real, double &scale_real_to_virtual,
                    double &offset_virtual_to_real, double &offset_real_to_virtual)
{
    double rg_virtual = max_virtual - min_virtual;
    double rg_viva = max_viva - min_viva;

    scale_virtual_to_real = rg_viva / rg_virtual;
    scale_real_to_virtual = rg_virtual / rg_viva;

    offset_virtual_to_real = min_viva - (min_virtual * scale_virtual_to_real);
    offset_real_to_virtual = min_virtual - (min_viva * scale_real_to_virtual);
}

void select_dt_neuron_model(const std::vector<double> &dts,
                            const std::vector<double> &pts,
                            double pts_live, double &dt, double &pts_burst)
{
    double aux = pts_live;
    double factor = 1.0;
    double intpart, fractpart;
    bool flag = false;

    dt = -1.0;
    pts_burst = -1.0;

    while (aux < pts[0])
    {
        aux = pts_live * factor;
        factor += 1.0;

        for (int i = pts.size() - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {
                dt = dts[i];
                pts_burst = pts[i];

                fractpart = std::modf(pts_burst / pts_live, &intpart);

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
        for (int i = pts.size() - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {
                dt = dts[i];
                pts_burst = pts[i];
                break;
            }
        }
    }
}

ScaledSignalResult scale_signal(
    const std::string &csv_path,
    int column_index,
    double csv_step,
    double start_time,
    double use_time,
    double observation_time,
    double freq,
    int integrator,
    bool check_drift,
    double sec_per_burst)
{

    ScaledSignalResult result;
    result.success = false;
    result.dt = -1.0;

    // Read CSV data
    std::vector<double> signal = read_csv_column(csv_path, column_index, start_time, use_time, csv_step);

    if (signal.empty())
    {
        return result;
    }

    // Calculate observation range
    size_t obs_points = static_cast<size_t>(observation_time / csv_step);
    if (obs_points > signal.size())
        obs_points = signal.size();

    std::vector<double> obs_signal(signal.begin(), signal.begin() + obs_points);

    // Determine min/max from observation
    double max_abs_real = -DBL_MAX;
    double min_abs_real = DBL_MAX;

    for (double val : obs_signal)
    {
        if (val > max_abs_real)
            max_abs_real = val;
        if (val < min_abs_real)
            min_abs_real = val;
    }

    double range = max_abs_real - min_abs_real;
    double percentage_min = 0.10;
    double percentage_max = 0.90;
    double min_rel_real = percentage_min * range + min_abs_real;
    double max_rel_real = percentage_max * range + min_abs_real;

    // Calculate signal period
    double external_firing_rate = signal_period(observation_time, obs_signal, max_rel_real, min_rel_real);

    if (sec_per_burst != -1.0)
    {
        external_firing_rate = sec_per_burst;
    }

    // Calculate dt and pts_burst
    double external_pts_per_burst = freq * external_firing_rate;
    double pts_burst = -1.0;

    if (integrator == RK4)
    {
        select_dt_neuron_model(HR_DTS_RK4, HR_PTS_RK4, external_pts_per_burst, result.dt, pts_burst);
    }

    if (result.dt == -1.0)
    {
        return result;
    }

    // Calculate s_points
    int s_points = static_cast<int>(pts_burst / external_pts_per_burst);
    if (s_points == 0)
        s_points = 1;

    // Initial scaling factors
    double min_abs_model = HINDMARSH_ROSE_MIN;
    double max_abs_model = HINDMARSH_ROSE_MAX;

    double scale_virtual_to_real, scale_real_to_virtual;
    double offset_virtual_to_real, offset_real_to_virtual;

    calcula_escala(min_abs_model, max_abs_model, min_abs_real, max_abs_real,
                   scale_virtual_to_real, scale_real_to_virtual,
                   offset_virtual_to_real, offset_real_to_virtual);

    // Apply vertical scaling to entire signal
    std::vector<double> scaled_signal;

    if (check_drift)
    {
        // Drift checking logic
        const int drift_n_burst = 2;
        int drift_counter = 0;
        double max_window = -999999.0;
        double min_window = 999999.0;
        double drift_aux_range = max_abs_real - min_abs_real;

        for (const auto& val : signal)
        {
            // Update window
            if ((min_window > val) && (val > (min_abs_real - drift_aux_range)))
            {
                min_window = val;
            }
            if ((max_window < val) && (val < (max_abs_real + drift_aux_range)))
            {
                max_window = val;
            }

            // Recalculate every drift_n_burst bursts
            if (drift_counter >= (drift_n_burst * external_pts_per_burst) &&
                max_window != -999999.0 && min_window != 999999.0)
            {

                calcula_escala(min_abs_model, max_abs_model, min_window, max_window,
                               scale_virtual_to_real, scale_real_to_virtual,
                               offset_virtual_to_real, offset_real_to_virtual);

                drift_aux_range = max_window - min_window;

                double per_min = 0.1, per_max = 0.1;
                if (min_window > 0)
                {
                    min_rel_real = min_window + (min_window * per_min);
                }
                else
                {
                    min_rel_real = min_window - (min_window * per_min);
                }

                if (max_window > 0)
                {
                    max_rel_real = max_window - (max_window * per_max);
                }
                else
                {
                    max_rel_real = max_window + (max_window * per_max);
                }

                max_window = -999999.0;
                min_window = 999999.0;
                drift_counter = 0;
            }

            drift_counter++;

            // Scale current point (vertical scaling)
            double scaled_val = val * scale_real_to_virtual + offset_real_to_virtual;
            scaled_signal.push_back(scaled_val);
        }
    }
    else
    {
        // No drift checking - simple scaling
        for (double val : signal)
        {
            double scaled_val = val * scale_real_to_virtual + offset_real_to_virtual;
            scaled_signal.push_back(scaled_val);
        }
    }

    // Perform linear interpolation (horizontal scaling)
    std::vector<double> interpolated_signal;

    for (size_t i = 0; i < scaled_signal.size() - 1; i++)
    {
        interpolated_signal.push_back(scaled_signal[i]);

        // Add s_points - 1 intermediate points
        for (int j = 1; j < s_points; j++)
        {
            double alpha = static_cast<double>(j) / s_points;
            double interp_val = scaled_signal[i] + alpha * (scaled_signal[i + 1] - scaled_signal[i]);
            interpolated_signal.push_back(interp_val);
        }
    }

    // Add last point
    interpolated_signal.push_back(scaled_signal.back());

    result.scaled_signal = interpolated_signal;
    result.success = true;

    return result;
}
