#include "scaling.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "utils.hpp"

using namespace kfr;

// Signal preprocessing for offline BO.
// Converts a raw external trace into:
// - scaled signal in model voltage domain
// - interpolation points for sub-step integration
// - selected dt / points_factor consistent with model burst dynamics

namespace SigPrivateConfig
{
    // Allowed relative mismatch between chosen points-per-burst and integer multiple.
    static constexpr double DT_SELECTION_TOLERANCE = 0.1;

    // Drift compensation percentages for rolling min/max windows.
    static constexpr double DRIFT_PERCENTAGE_MIN = 0.1;
    static constexpr double DRIFT_PERCENTAGE_MAX = 0.1;
    // Number of detected bursts between drift recalibrations.
    static constexpr size_t DRIFT_N_BURST = 2;
}

struct ScalingFactors
{
    // Affine transform: virtual = real * scale + offset.
    double scale_real_to_virtual;
    double offset_real_to_virtual;
};

struct DTSelection
{
    // Selected integration step and matched points-per-burst table value.
    double dt;
    double pts_burst;
};

struct SigStats
{
    // Absolute/relative stats and estimated period from observation window.
    double real_abs_min;
    double real_abs_max;
    double real_rel_min;
    double real_rel_max;
    double sig_period;
};

static inline void read_csv_column(univector<double> &data, const std::string &csv_path, size_t column_idx,
                                   size_t start_idx, size_t num_points)
{
    // Read one CSV column in [start_idx, start_idx + num_points) into `data`.
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
            // Reuse stringstream buffer line-by-line to avoid extra allocations.
            ss.str(line);
            ss.clear();
            size_t current_col = 0;

            while (std::getline(ss, val, ','))
            {
                if (current_col == column_idx)
                {
                    // Parse only requested column and append it.
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
    // Estimate signal period by counting low->high crossings with hysteresis.
    // `th_up` and `th_on` reduce false transitions in noisy traces.
    bool up = (sig.front() > th_up);
    // Number of completed upward transitions in the observation window.
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

    // period = total_time / num_cycles.
    return 1.0 / (changes / observation_time);
}

static ScalingFactors calculate_scaling(double virtual_min, double virtual_max,
                                        double real_min, double real_max)
{
    // Compute affine mapping from real-domain values to model-domain values.
    // Endpoint mapping:
    // real_min -> virtual_min, real_max -> virtual_max.
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
    // Select dt from lookup tables so model points/burst approximates signal points/burst
    // with near-integer interpolation factor.
    // `aux` is the target points-per-burst we try to match in lookup table space.
    double aux = pts_real;
    // `factor` explores integer multiples of pts_real to find near-integer ratios.
    double factor = 1.0;
    double intpart, fractpart;
    bool flag = false;

    constexpr double INVALID_DT = SigConstants::INVALID_DT;

    double dt_candidate = INVALID_DT;
    double pts_burst_candidate = SigConstants::INVALID_PTS;

    // First pass: search for a candidate whose ratio pts_model/pts_real is near-integer.
    while (aux < pts[0])
    {
        // Try next multiple of the real signal points-per-burst.
        aux = pts_real * factor;
        factor += 1.0;

        // Iterate from largest points-per-burst to smallest.
        for (size_t i = N - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {
                // Keep the closest table entry above current target.
                dt_candidate = dts[i];
                pts_burst_candidate = pts[i];

                fractpart = std::modf(pts_burst_candidate / pts_real, &intpart);

                // Prefer candidates with small fractional mismatch.
                if (fractpart <= SigPrivateConfig::DT_SELECTION_TOLERANCE * intpart)
                {
                    flag = true;
                }
                // Always break before i underflows (size_t loop).
                break;
            }
        }

        if (flag)
        {
            break;
        }
    }

    // Fallback: if near-integer criterion fails, keep best feasible candidate.
    if (!flag)
    {
        for (size_t i = N - 1; i >= 0; i--)
        {
            if (pts[i] > aux)
            {
                // Best feasible fallback even if ratio is not near-integer.
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
    // Dispatch dt/points table by neuron model and numeric integrator.
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
    // Re-estimate scaling factors from a local window to compensate baseline drift.
    ScalingFactors factors = calculate_scaling(model_abs_min, model_abs_max, window_min, window_max);

    constexpr double DRIFT_PERCENTAGE_MIN = SigPrivateConfig::DRIFT_PERCENTAGE_MIN;
    constexpr double DRIFT_PERCENTAGE_MAX = SigPrivateConfig::DRIFT_PERCENTAGE_MAX;

    // Update relative thresholds with sign-aware margins.
    if (window_min > 0)
    {
        // If positive, push threshold slightly upward.
        stats.real_rel_min = window_min + (window_min * DRIFT_PERCENTAGE_MIN);
    }
    else
    {
        // If negative, move threshold to more negative value.
        stats.real_rel_min = window_min - (window_min * DRIFT_PERCENTAGE_MIN);
    }

    if (window_max > 0)
    {
        // If positive, pull threshold slightly downward.
        stats.real_rel_max = window_max - (window_max * DRIFT_PERCENTAGE_MAX);
    }
    else
    {
        // If negative, move threshold to less negative value.
        stats.real_rel_max = window_max + (window_max * DRIFT_PERCENTAGE_MAX);
    }

    return factors;
}

static inline SigStats init_stats(const univector<double> &sig, size_t obs_points, double csv_step)
{
    // Extract observation statistics used for period estimation and scaling.
    const double observation_time_to_use = obs_points * csv_step;

    SigStats stats;

    // Observation window is always taken from the beginning of the loaded segment.
    const univector_ref<const double> obs_sig = sig.slice(0, obs_points);
    const double abs_max = maxof(obs_sig);
    const double abs_min = minof(obs_sig);

    stats.real_abs_min = abs_min;
    stats.real_abs_max = abs_max;

    const double range = abs_max - abs_min;
    // Relative thresholds (10%-90% of range) for robust transition detection.
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
    // Scaling pipeline:
    // 1) Read CSV segment.
    // 2) Estimate period and choose dt/points factor from calibration tables.
    // 3) Rescale to model voltage range (with optional drift compensation).
    // 4) Build per-interval interpolated points for inner integration steps.
    if (csv_step <= 0 || use_time <= 0 || observation_time <= 0 || start_time < 0 || column_idx < 0 || csv_path.empty())
    {
        throw std::invalid_argument("Invalid arguments: csv_step, use_time, observation_time must be positive, start_time and column_idx non-negative, csv_path non-empty");
    }

    ScaledSigResult result;

    univector<double> &sig = result.sig;
    univector<double> &interpolated_points = result.interpolated_points;

    // Time-to-index conversion (floor behavior from static_cast).
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
    // Need enough points for both:
    // - observation window (period estimation)
    // - use window (actual BO signal)
    const size_t read_points = std::max(use_points, obs_points);

    read_csv_column(sig, csv_path, column_idx, start_idx, read_points);

    // read_csv_column may return fewer points if file is shorter than requested.
    size_t sig_size = sig.size();
    if (sig_size == 0)
    {
        throw std::runtime_error("No data read from CSV file");
    }

    SigStats stats = init_stats(sig, obs_points, csv_step);

    // Keep only the segment requested for BO usage.
    if (sig_size > use_points)
    {
        sig.resize(use_points);
        sig_size = use_points;
    }

    // External points per burst at original sampling step.
    const double external_pts_per_burst = stats.sig_period / csv_step;

    std::optional<DTSelection> selection;
    double model_abs_min, model_abs_max;
    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        // Target model voltage interval for affine rescaling.
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
        // Could not find a valid dt candidate for this signal/model pair.
        return std::nullopt;
    }

    result.dt = selection->dt;

    // Integration sub-steps per original sample.
    size_t s_points = static_cast<size_t>(selection->pts_burst / external_pts_per_burst);
    // Keep at least one step per sample in degenerate cases.
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
        // Online-style drift correction: periodically update scaling with local windows.
        size_t drift_counter = 0;
        constexpr double DOUBLE_MAX = GeneralConstants::DOUBLE_MAX;
        constexpr double DOUBLE_MIN = GeneralConstants::DOUBLE_MIN;
        double window_max = DOUBLE_MIN;
        double window_min = DOUBLE_MAX;
        const double drift_aux_range = real_abs_max - real_abs_min;

        for (double &val : sig)
        {
            // Keep local extrema within a bounded neighborhood to reject outliers.
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
                // Recalibrate every N bursts (in source-sample units).
                drift_counter = 0;

                // Refit affine transform periodically from local window extrema.
                factors = fix_drift(model_abs_min, model_abs_max, window_min, window_max, stats);
                scale_real_to_virtual = factors.scale_real_to_virtual;
                offset_real_to_virtual = factors.offset_real_to_virtual;

                window_max = DOUBLE_MIN;
                window_min = DOUBLE_MAX;
            }

            drift_counter++;

            // Apply current affine mapping.
            val = val * scale_real_to_virtual + offset_real_to_virtual;
        }
    }
    else
    {
        // Global affine scaling for the whole signal.
        sig = (sig * scale_real_to_virtual) + offset_real_to_virtual;
    }

    // Precompute interpolated values between consecutive samples.
    const size_t interpolated_size = (sig_size - 1) * (s_points - 1);
    interpolated_points.reserve(interpolated_size);
    const double *sig_ptr = sig.data();

    for (size_t i = 0; i < sig_size - 1; i++)
    {
        // Create (s_points - 1) inner linear samples per original interval.
        for (double j = 1.0; j < s_points; j++)
        {
            // j is double to ensure floating-point interpolation ratio.
            double alpha = j / s_points;
            double interp_val = sig_ptr[i] + (alpha * (sig_ptr[i + 1] - sig_ptr[i]));
            interpolated_points.push_back(interp_val);
        }
    }

    return result;
}
