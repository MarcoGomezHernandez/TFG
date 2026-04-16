#include "scaling.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "utils.hpp"

using namespace kfr;

namespace SigPrivateConfig
{
    static constexpr double DT_SELECTION_TOLERANCE = 0.1;

    static constexpr double DRIFT_PERCENTAGE_MIN = 0.1;
    static constexpr double DRIFT_PERCENTAGE_MAX = 0.1;
    static constexpr size_t DRIFT_N_BURST = 2;
}

struct ScalingFactors
{
    double scale_real_to_virtual;
    double offset_real_to_virtual;
};

struct DTSelection
{
    double dt;
    double pts_burst;
};

struct SigStats
{
    double real_abs_min;
    double real_abs_max;
    double real_rel_min;
    double real_rel_max;
    double sig_period;
};

static inline void read_csv_column(univector<double> &data, const std::string &csv_path, size_t column_idx,
                                   size_t start_idx, size_t num_points)
{
    std::ifstream file(csv_path);

    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open CSV file: " + csv_path);
    }

    data.reserve(num_points);

    const size_t end_idx = start_idx + num_points;

    std::string line;
    std::string val;
    std::stringstream ss;
    size_t current_line = 0;

    while (current_line < end_idx && std::getline(file, line))
    {
        if (current_line >= start_idx)
        {
            ss.str(line);
            ss.clear();
            size_t current_col = 0;

            while (std::getline(ss, val, ','))
            {
                if (current_col == column_idx)
                {
                    data.push_back(std::stod(val));
                    break;
                }
                current_col++;
            }
        }

        current_line++;
    }

    const size_t data_size = data.size();

    if (data_size < num_points)
    {
        data.shrink_to_fit();
        std::cout << "Warning: Fewer data points read (" << data_size << ") than expected (" << num_points << ") from CSV file: " << csv_path << std::endl;
    }

    file.close();
}

static inline double sig_period(double observation_time, const univector_ref<const double> &sig,
                                double th_up, double th_on)
{
    bool up = (sig.front() > th_up);
    double changes = 0.0;

    for (const double val : sig)
    {
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

    return 1.0 / (changes / observation_time);
}

static ScalingFactors calculate_scaling(double virtual_min, double virtual_max,
                                        double real_min, double real_max)
{
    const double virtual_range = virtual_max - virtual_min;
    const double real_range = real_max - real_min;

    ScalingFactors factors;
    const double scale_real_to_virtual = virtual_range / real_range;
    factors.scale_real_to_virtual = scale_real_to_virtual;
    factors.offset_real_to_virtual = virtual_min - (real_min * scale_real_to_virtual);

    return factors;
}

template <size_t N>
static inline std::optional<DTSelection> select_dt_neuron_model(const std::array<double, N> &dts,
                                                                const std::array<double, N> &pts,
                                                                double pts_real)
{
    double aux = pts_real;
    double factor = 1.0;
    double intpart, fractpart;
    bool flag = false;

    constexpr double INVALID_DT = SigConstants::INVALID_DT;

    double dt_candidate = INVALID_DT;
    double pts_burst_candidate = SigConstants::INVALID_PTS;

    while (aux < pts[0])
    {
        aux = pts_real * factor;
        factor += 1.0;

        for (size_t i = N - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {
                dt_candidate = dts[i];
                pts_burst_candidate = pts[i];

                fractpart = std::modf(pts_burst_candidate / pts_real, &intpart);

                if (fractpart <= SigPrivateConfig::DT_SELECTION_TOLERANCE * intpart)
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

    if (dt_candidate == INVALID_DT)
    {
        return std::nullopt;
    }

    DTSelection selection;
    selection.dt = dt_candidate;
    selection.pts_burst = pts_burst_candidate;

    return selection;
}

static inline std::optional<DTSelection> set_pts_burst(NeuronModel model, NumericIntegrator integrator, double pts_real)
{
    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        if (integrator == NumericIntegrator::RK4)
        {
            return select_dt_neuron_model(HindmarshRose::DTS_RK4, HindmarshRose::PTS_RK4, pts_real);
        }
        else
        {
            throw std::invalid_argument("Unsupported integrator");
        }
    }
    else
    {
        throw std::invalid_argument("Unsupported neuron model");
    }
}

static inline ScalingFactors fix_drift(double model_abs_min, double model_abs_max, double window_min, double window_max, SigStats &stats)
{
    ScalingFactors factors = calculate_scaling(model_abs_min, model_abs_max, window_min, window_max);

    constexpr double DRIFT_PERCENTAGE_MIN = SigPrivateConfig::DRIFT_PERCENTAGE_MIN;
    constexpr double DRIFT_PERCENTAGE_MAX = SigPrivateConfig::DRIFT_PERCENTAGE_MAX;

    if (window_min > 0)
    {
        stats.real_rel_min = window_min + (window_min * DRIFT_PERCENTAGE_MIN);
    }
    else
    {
        stats.real_rel_min = window_min - (window_min * DRIFT_PERCENTAGE_MIN);
    }

    if (window_max > 0)
    {
        stats.real_rel_max = window_max - (window_max * DRIFT_PERCENTAGE_MAX);
    }
    else
    {
        stats.real_rel_max = window_max + (window_max * DRIFT_PERCENTAGE_MAX);
    }

    return factors;
}

static inline SigStats init_stats(const univector<double> &sig, size_t obs_points, double csv_step)
{
    const double observation_time_to_use = obs_points * csv_step;

    SigStats stats;

    const univector_ref<const double> obs_sig = sig.slice(0, obs_points);
    const double abs_max = maxof(obs_sig);
    const double abs_min = minof(obs_sig);

    stats.real_abs_min = abs_min;
    stats.real_abs_max = abs_max;

    const double range = abs_max - abs_min;
    const double real_rel_min = SigPublicConfig::SIG_PERCENTAGE_MIN * range + abs_min;
    const double real_rel_max = SigPublicConfig::SIG_PERCENTAGE_MAX * range + abs_min;
    stats.real_rel_min = real_rel_min;
    stats.real_rel_max = real_rel_max;

    stats.sig_period = sig_period(observation_time_to_use, obs_sig, real_rel_max, real_rel_min);

    return stats;
}

std::optional<ScaledSigResult> scale_sig(
    const std::string &csv_path,
    size_t column_idx,
    double csv_step,
    double start_time,
    double use_time,
    double observation_time,
    NumericIntegrator integrator,
    NeuronModel model,
    bool check_drift)
{
    if (csv_step <= 0 || use_time <= 0 || observation_time <= 0 || start_time < 0 || column_idx < 0 || csv_path.empty())
    {
        throw std::invalid_argument("Invalid arguments: csv_step, use_time, observation_time must be positive, start_time and column_idx non-negative, csv_path non-empty");
    }

    ScaledSigResult result;

    univector<double> &sig = result.sig;
    univector<double> &interpolated_points = result.interpolated_points;

    const size_t start_idx = static_cast<size_t>(start_time / csv_step);
    const size_t use_points = static_cast<size_t>(use_time / csv_step);
    if (use_points == 0)
    {
        throw std::invalid_argument("use_time is too short to read any points with given csv_step");
    }

    const size_t obs_points = static_cast<size_t>(observation_time / csv_step);
    if (obs_points == 0)
    {
        throw std::invalid_argument("observation_time is too short to read any points with given csv_step");
    }
    const size_t read_points = std::max(use_points, obs_points);

    read_csv_column(sig, csv_path, column_idx, start_idx, read_points);

    size_t sig_size = sig.size();
    if (sig_size == 0)
    {
        throw std::runtime_error("No data read from CSV file");
    }

    SigStats stats = init_stats(sig, obs_points, csv_step);

    if (sig_size > use_points)
    {
        sig.resize(use_points);
        sig_size = use_points;
    }

    const double external_pts_per_burst = stats.sig_period / csv_step;

    std::optional<DTSelection> selection;
    double model_abs_min, model_abs_max;
    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        model_abs_min = HindmarshRose::MIN;
        model_abs_max = HindmarshRose::MAX;
        selection = set_pts_burst(model, integrator, external_pts_per_burst);
    }
    else
    {
        throw std::invalid_argument("Unsupported neuron model");
    }

    if (!selection)
    {
        return std::nullopt;
    }

    result.dt = selection->dt;

    size_t s_points = static_cast<size_t>(selection->pts_burst / external_pts_per_burst);
    if (s_points == 0)
        s_points = 1;

    result.points_factor = s_points;

    double &real_abs_min = stats.real_abs_min;
    double &real_abs_max = stats.real_abs_max;

    ScalingFactors factors = calculate_scaling(model_abs_min, model_abs_max, real_abs_min, real_abs_max);
    double scale_real_to_virtual = factors.scale_real_to_virtual;
    double offset_real_to_virtual = factors.offset_real_to_virtual;

    if (check_drift)
    {
        size_t drift_counter = 0;
        constexpr double DOUBLE_MAX = GeneralConstants::DOUBLE_MAX;
        constexpr double DOUBLE_MIN = GeneralConstants::DOUBLE_MIN;
        double window_max = DOUBLE_MIN;
        double window_min = DOUBLE_MAX;
        const double drift_aux_range = real_abs_max - real_abs_min;

        for (double &val : sig)
        {
            if ((window_min > val) && (val > (real_abs_min - drift_aux_range)))
            {
                window_min = val;
            }
            if ((window_max < val) && (val < (real_abs_max + drift_aux_range)))
            {
                window_max = val;
            }

            if (drift_counter >= (SigPrivateConfig::DRIFT_N_BURST * external_pts_per_burst) &&
                window_max != DOUBLE_MIN && window_min != DOUBLE_MAX)
            {
                drift_counter = 0;

                factors = fix_drift(model_abs_min, model_abs_max, window_min, window_max, stats);
                scale_real_to_virtual = factors.scale_real_to_virtual;
                offset_real_to_virtual = factors.offset_real_to_virtual;

                window_max = DOUBLE_MIN;
                window_min = DOUBLE_MAX;
            }

            drift_counter++;

            val = val * scale_real_to_virtual + offset_real_to_virtual;
        }
    }
    else
    {
        sig = (sig * scale_real_to_virtual) + offset_real_to_virtual;
    }

    const size_t interpolated_size = (sig_size - 1) * (s_points - 1);
    interpolated_points.reserve(interpolated_size);
    const double *sig_ptr = sig.data();

    for (size_t i = 0; i < sig_size - 1; i++)
    {
        for (double j = 1.0; j < s_points; j++)
        {
            double alpha = j / s_points;
            double interp_val = sig_ptr[i] + (alpha * (sig_ptr[i + 1] - sig_ptr[i]));
            interpolated_points.push_back(interp_val);
        }
    }

    return result;
}
