#include "scaling.hpp"
#include <kfr/all.hpp>
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
    bool success;
};

struct SigStats
{
    double min_abs_real;
    double max_abs_real;
    double min_rel_real;
    double max_rel_real;
    double period_sig;
};

void read_csv_column(univector<double> &data, const std::string &csv_path, size_t column_i,
                     size_t start_i, size_t num_points)
{
    std::ifstream file(csv_path);

    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open CSV file: " + csv_path);
    }

    data.reserve(num_points);

    const size_t end_i = start_i + num_points;

    std::string line;
    std::string val;
    std::stringstream ss;
    size_t current_line = 0;

    while (current_line < end_i && std::getline(file, line))
    {
        if (current_line >= start_i)
        {
            ss.str(line);
            ss.clear();
            size_t current_col = 0;

            while (std::getline(ss, val, ','))
            {
                if (current_col == column_i)
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

double sig_period(double tiempo_observacion, const univector_ref<double> &sig,
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

    return 1.0 / (changes / tiempo_observacion);
}

ScalingFactors calcula_escala(double min_virtual, double max_virtual,
                              double min_viva, double max_viva)
{
    const double rg_virtual = max_virtual - min_virtual;
    const double rg_viva = max_viva - min_viva;

    ScalingFactors factors;
    const double scale_real_to_virtual = rg_virtual / rg_viva;
    factors.scale_real_to_virtual = scale_real_to_virtual;
    factors.offset_real_to_virtual = min_virtual - (min_viva * scale_real_to_virtual);

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

    constexpr double INVALID_DT = SigConstants::INVALID_DT;

    double dt_candidate = INVALID_DT;
    double pts_burst_candidate = SigConstants::INVALID_PTS;

    while (aux < pts[0])
    {
        aux = pts_live * factor;
        factor += 1.0;

        for (size_t i = N - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {
                dt_candidate = dts[i];
                pts_burst_candidate = pts[i];

                fractpart = std::modf(pts_burst_candidate / pts_live, &intpart);

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

    DTSelection selection;

    selection.success = (dt_candidate != INVALID_DT);
    selection.dt = dt_candidate;
    selection.pts_burst = pts_burst_candidate;

    return selection;
}

inline DTSelection set_pts_burst(NeuronModel model, NumericIntegrator integrator, double pts_live)
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

inline ScalingFactors fix_drift(double min_abs_model, double max_abs_model, double min_window, double max_window, SigStats &stats)
{
    ScalingFactors factors = calcula_escala(min_abs_model, max_abs_model, min_window, max_window);

    constexpr double DRIFT_PERCENTAGE_MIN = SigPrivateConfig::DRIFT_PERCENTAGE_MIN;
    constexpr double DRIFT_PERCENTAGE_MAX = SigPrivateConfig::DRIFT_PERCENTAGE_MAX;

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

SigStats ini_recibido(const univector<double> &sig, size_t obs_points, double csv_step)
{
    const double observation_time_to_use = obs_points * csv_step;

    SigStats stats;

    const univector_ref<double> obs_sig = sig.slice(0, obs_points);
    const double max_abs = maxof(obs_sig);
    const double min_abs = minof(obs_sig);

    stats.min_abs_real = min_abs;
    stats.max_abs_real = max_abs;

    const double range = max_abs - min_abs;
    const double min_rel_real = SigPublicConfig::SIG_PERCENTAGE_MIN * range + min_abs;
    const double max_rel_real = SigPublicConfig::SIG_PERCENTAGE_MAX * range + min_abs;
    stats.min_rel_real = min_rel_real;
    stats.max_rel_real = max_rel_real;

    stats.period_sig = sig_period(observation_time_to_use, obs_sig, max_rel_real, min_rel_real);

    return stats;
}

ScaledSigResult scale_sig(
    const std::string &csv_path,
    size_t column_i,
    double csv_step,
    double start_time,
    double use_time,
    double observation_time,
    NumericIntegrator integrator,
    NeuronModel model,
    bool check_drift)
{
    if (csv_step <= 0 || use_time <= 0 || observation_time <= 0 || start_time < 0 || column_i < 0 || csv_path.empty())
    {
        throw std::runtime_error("Invalid arguments: csv_step, use_time, observation_time must be positive, start_time and column_i non-negative, csv_path non-empty");
    }

    ScaledSigResult result;

    univector<double> &sig = result.sig;
    univector<double> &interpolated_points = result.interpolated_points;

    const size_t start_i = static_cast<size_t>(start_time / csv_step);
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

    read_csv_column(sig, csv_path, column_i, start_i, read_points);

    size_t sig_size = sig.size();
    if (sig_size == 0)
    {
        throw std::runtime_error("No data read from CSV file");
    }

    SigStats stats = ini_recibido(sig, obs_points, csv_step);

    if (sig_size > use_points)
    {
        sig.resize(use_points);
        sig_size = use_points;
    }

    const double external_pts_per_burst = stats.period_sig / csv_step;

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

    if (!selection.success)
    {
        result.success = false;
        sig.clear();
        result.dt = SigConstants::INVALID_DT;
        return result;
    }

    result.dt = selection.dt;
    result.pts_burst_real = external_pts_per_burst;

    size_t s_points = static_cast<size_t>(selection.pts_burst / external_pts_per_burst);
    if (s_points == 0)
        s_points = 1;

    result.points_factor = s_points;

    double &min_abs_real = stats.min_abs_real;
    double &max_abs_real = stats.max_abs_real;

    ScalingFactors factors = calcula_escala(min_abs_model, max_abs_model, min_abs_real, max_abs_real);
    double scale_real_to_virtual = factors.scale_real_to_virtual;
    double offset_real_to_virtual = factors.offset_real_to_virtual;

    if (check_drift)
    {
        size_t drift_counter = 0;
        constexpr double DOUBLE_MAX = GeneralConstants::DOUBLE_MAX;
        constexpr double DOUBLE_MIN = GeneralConstants::DOUBLE_MIN;
        double max_window = DOUBLE_MIN;
        double min_window = DOUBLE_MAX;
        const double drift_aux_range = max_abs_real - min_abs_real;

        for (double &val : sig)
        {
            if ((min_window > val) && (val > (min_abs_real - drift_aux_range)))
            {
                min_window = val;
            }
            if ((max_window < val) && (val < (max_abs_real + drift_aux_range)))
            {
                max_window = val;
            }

            if (drift_counter >= (SigPrivateConfig::DRIFT_N_BURST * external_pts_per_burst) &&
                max_window != DOUBLE_MIN && min_window != DOUBLE_MAX)
            {
                drift_counter = 0;

                factors = fix_drift(min_abs_model, max_abs_model, min_window, max_window, stats);
                scale_real_to_virtual = factors.scale_real_to_virtual;
                offset_real_to_virtual = factors.offset_real_to_virtual;

                max_window = DOUBLE_MIN;
                min_window = DOUBLE_MAX;
            }

            drift_counter++;

            val = val * scale_real_to_virtual + offset_real_to_virtual;
        }
    }
    else
    {
        sig = sig * scale_real_to_virtual + offset_real_to_virtual;
    }

    const size_t interpolated_size = (sig_size - 1) * (s_points - 1);
    interpolated_points.reserve(interpolated_size);

    for (size_t i = 0; i < sig_size - 1; i++)
    {
        for (double j = 1.0; j < s_points; j++)
        {
            double alpha = j / s_points;
            double interp_val = sig[i] + (alpha * (sig[i + 1] - sig[i]));
            interpolated_points.push_back(interp_val);
        }
    }

    result.success = true;

    return result;
}
