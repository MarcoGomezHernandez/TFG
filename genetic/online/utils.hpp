#ifndef EVALUATION_UTILS_H
#define EVALUATION_UTILS_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <Eigen/Core>
#include <kfr/all.hpp>

namespace BOPublicConfig
{
    // Extra scaling applied to the allowed distance between expected current bounds.
    inline constexpr double EXPECTED_I_MARGIN_RANGE_FACTOR = 0.8;

    // Active-wait sleep duration used by RT/NRT handshakes (milliseconds).
    inline constexpr double ACTIVE_WAIT_MS = 10.0;

    // Margins used when expanding voltage bounds for BO parameter initialization.
    inline constexpr double V_MARGIN_FACTOR = 0.125;
    // Minimum position of V_fast inside the observed pre-synaptic voltage range.
    inline constexpr double V_FAST_MIN_FACTOR = 0.25;

    // Margin around the user-provided expected current range (min/max) used for scoring.
    inline constexpr double EXPECTED_I_MARGIN_FACTOR = 0.5;

    // Sigmoid slope factors for fast/slow components (scaled by 1 / V_range).
    inline constexpr double S_FAST_MIN_FACTOR = 4.12;
    inline constexpr double S_FAST_MAX_FACTOR = 13.72;
    inline constexpr double S_SLOW_MIN_FACTOR = 1.72;
    inline constexpr double S_SLOW_MAX_FACTOR = 3.43;

    // k1 is scaled by cutoff frequency (fc) to keep time constants consistent.
    inline constexpr double K1_FACTOR = 3.33;
    // Lower bound factor for k1 (prevents degenerate near-zero dynamics).
    inline constexpr double K1_MIN_FACTOR = 0.000001;

    // How far/near the synaptic reversal potential (E_syn) is allowed from V_post.
    // (Scaled by the observed post-synaptic voltage range.)
    inline constexpr double E_SYN_FAR_TERM = 3.86;
    inline constexpr double E_SYN_NEAR_TERM = 0.2;

    // Lower bound multiplier for conductances (relative to computed max).
    inline constexpr double G_MIN_FACTOR = 0.001;

    // Bounds for R = k2/k1 (dimensionless ratio in the slow gating ODE).
    inline constexpr double R_MAX = 100.0;
    inline constexpr double R_MIN = 0.01;
}

namespace BOPublicConstants
{
    // Small positive values used to avoid log(0) and division by ~0.
    inline constexpr double SMALL_LOG = std::numeric_limits<double>::min();
    inline constexpr double SMALL_DIVISOR = std::numeric_limits<double>::epsilon();
    inline constexpr double NEGATIVE_SMALL_DIVISOR = -SMALL_DIVISOR;
}

// Exception used to abort a Limbo evaluation/optimization cleanly.
struct StopEvaluation
{
};

// Standard logistic activation used by the chemical synapse model.
inline double chemical_sigmoid(double s,
                               double v_threshold,
                               double v_pre)
{
    return 1.0 / (1.0 + std::exp(s * (v_threshold - v_pre)));
}

// Computes a normalized "max distance" used when comparing expected vs observed min/max.
// This includes margins around the expected current range.
inline double calculate_expected_i_max_dist(double expected_i_min,
                                     double expected_i_max)
{
    const double range = expected_i_max - expected_i_min;
    return (range + (range * BOPublicConfig::EXPECTED_I_MARGIN_FACTOR * 2.0)) *
           BOPublicConfig::EXPECTED_I_MARGIN_RANGE_FACTOR;
}

inline double safe_divisor(double divisor)
{
    // Prevent divisions from blowing up when divisor is ~0.
    return std::abs(divisor) < BOPublicConstants::SMALL_DIVISOR
               ? (divisor < 0.0 ? BOPublicConstants::NEGATIVE_SMALL_DIVISOR : BOPublicConstants::SMALL_DIVISOR)
               : divisor;
}

// Parameter decoding ranges for Limbo BO.
// Each ParamRange maps x\in[0,1] to a physical parameter via: value = min + x*range.
struct BOParamRanges
{
    struct ParamRange
    {
        // Lower bound and (max-min) span for decoding.
        double min;
        double range;

        ParamRange() = default;

        ParamRange(double min, double max)
            : min(min), range(max - min)
        {
        }
    };

    // Decoding ranges (some are log-domain: log_k1/log_R/log_g_*).
    ParamRange s_fast;
    ParamRange s_slow;
    ParamRange e_syn;
    ParamRange log_k1;
    ParamRange log_R;
    ParamRange log_g_fast;
    ParamRange log_g_slow;
    ParamRange v_fast;
    ParamRange v_slow;

    // Initialize all ranges for a single synapse direction.
    // Inputs:
    // - v_pre_* / v_post_*: observed voltage bounds used to constrain thresholds and E_syn
    // - expected_i_*: user target current bounds to match (for scoring and g bounds)
    // - use_i_fast/use_i_slow: enable/disable subsets of parameters
    // - search_phase: flips the allowed E_syn side relative to V_post
    // - fc: cutoff frequency used to scale slow dynamics (k1)
    void init(double v_pre_min,
              double v_pre_max,
              double v_post_min,
              double v_post_max,
              double expected_i_min,
              double expected_i_max,
              unsigned int use_i_fast,
              unsigned int use_i_slow,
              unsigned int search_phase,
              double fc)
    {
        constexpr double G_MIN_FACTOR = BOPublicConfig::G_MIN_FACTOR;
        constexpr double R_MIN = BOPublicConfig::R_MIN;
        constexpr double SMALL_LOG = BOPublicConstants::SMALL_LOG;
        constexpr double R_MAX = BOPublicConfig::R_MAX;

        // Pre/post voltage spans.
        const double v_pre_range = v_pre_max - v_pre_min;
        const double v_post_range = v_post_max - v_post_min;

        // E_syn bounds are positioned "near" or "far" from V_post depending on phase.
        const double e_syn_far_final_term = v_post_range * BOPublicConfig::E_SYN_FAR_TERM;
        const double e_syn_near_final_term = v_post_range * BOPublicConfig::E_SYN_NEAR_TERM;
        double e_syn_max = 0.0, e_syn_min = 0.0;
        if (search_phase)
        {
            e_syn_min = v_post_max + e_syn_near_final_term;
            e_syn_max = v_post_max + e_syn_far_final_term;
        }
        else
        {
            e_syn_min = v_post_min - e_syn_far_final_term;
            e_syn_max = v_post_min - e_syn_near_final_term;
        }
        e_syn = ParamRange(e_syn_min, e_syn_max);

        // Expand V_pre bounds a bit to avoid corner solutions and edge effects.
        const double v_pre_margin = v_pre_range * BOPublicConfig::V_MARGIN_FACTOR;
        const double v_pre_margin_min = v_pre_min - v_pre_margin;
        const double v_pre_margin_max = v_pre_max + v_pre_margin;

        // Expand expected current range to allow some slack during optimization.
        const double expected_i_margin = (expected_i_max - expected_i_min) * BOPublicConfig::EXPECTED_I_MARGIN_FACTOR;
        const double expected_i_margin_min = expected_i_min - expected_i_margin;
        const double expected_i_margin_max = expected_i_max + expected_i_margin;
        const double safe_v_pre_range = safe_divisor(v_pre_range);

        if (use_i_fast)
        {
            // Fast component: threshold constrained to upper part of V_pre, and steep sigmoid.
            v_fast = ParamRange(v_pre_min + (v_pre_range * BOPublicConfig::V_FAST_MIN_FACTOR),
                                v_pre_margin_max);
            const double s_fast_max = BOPublicConfig::S_FAST_MAX_FACTOR / safe_v_pre_range;
            s_fast = ParamRange(BOPublicConfig::S_FAST_MIN_FACTOR / safe_v_pre_range,
                                s_fast_max);

            // Bound g_fast by ensuring the expected current can be achieved at extremes.
            const double sigmoid_fast_max = chemical_sigmoid(s_fast_max, v_fast.min, v_pre_max);
            const double g_fast_max = std::max(std::abs(expected_i_margin_max / safe_divisor((v_post_max - e_syn_min) * sigmoid_fast_max)),
                                               std::abs(expected_i_margin_min / safe_divisor((v_post_min - e_syn_max) * sigmoid_fast_max)));
            const double g_fast_min = g_fast_max * G_MIN_FACTOR;
            log_g_fast = ParamRange(std::log(g_fast_min == 0.0 ? SMALL_LOG : g_fast_min),
                                    std::log(g_fast_max == 0.0 ? SMALL_LOG : g_fast_max));
        }

        if (use_i_slow)
        {
            // Slow component: k1 scales with fc; k2 = k1 * R with bounded ratio R.
            const double k1_max = BOPublicConfig::K1_FACTOR * fc;
            const double k1_min = k1_max * BOPublicConfig::K1_MIN_FACTOR;

            // Slow threshold can cover the whole (margined) V_pre range.
            v_slow = ParamRange(v_pre_margin_min,
                                v_pre_margin_max);
            const double s_slow_max = BOPublicConfig::S_SLOW_MAX_FACTOR / safe_v_pre_range;
            s_slow = ParamRange(BOPublicConfig::S_SLOW_MIN_FACTOR / safe_v_pre_range,
                                s_slow_max);

            // Log-domain ranges improve numerical conditioning during BO.
            log_k1 = ParamRange(std::log(k1_min == 0.0 ? SMALL_LOG : k1_min),
                                std::log(k1_max == 0.0 ? SMALL_LOG : k1_max));

            log_R = ParamRange(std::log(R_MIN == 0.0 ? SMALL_LOG : R_MIN),
                               std::log(R_MAX == 0.0 ? SMALL_LOG : R_MAX));

            // Conservative bound for the maximum achievable m (steady-state gate).
            const double k2_min = k1_min * R_MIN;

            const double m_max_term = k1_max * chemical_sigmoid(s_slow_max, v_slow.min, v_pre_max);
            const double m_max = m_max_term / safe_divisor(m_max_term + k2_min);

            // Bound g_slow similarly to g_fast, using m_max as the worst-case activation.
            const double g_slow_max = std::max(std::abs(expected_i_margin_max / safe_divisor((v_post_max - e_syn_min) * m_max)),
                                               std::abs(expected_i_margin_min / safe_divisor((v_post_min - e_syn_max) * m_max)));
            const double g_slow_min = g_slow_max * G_MIN_FACTOR;

            log_g_slow = ParamRange(std::log(g_slow_min == 0.0 ? SMALL_LOG : g_slow_min),
                                    std::log(g_slow_max == 0.0 ? SMALL_LOG : g_slow_max));
        }
    }

    BOParamRanges() = default;
};

struct EvaluationPadBuffers
{
    // Scratch buffers used to build padded segments for zero-phase filtering.
    kfr::univector<double> padded_buff_12;
    kfr::univector<double> padded_buff_21;

    EvaluationPadBuffers(size_t padded_buff_size_12,
                         size_t padded_buff_size_21)
        : padded_buff_12(padded_buff_size_12),
          padded_buff_21(padded_buff_size_21)
    {
    }
};

struct StopFunctor
{
    // Shared stop flag (owned by module) used as an external stop criterion in Limbo.
    static inline std::atomic<bool> *stop_BO_ptr = nullptr;

    template <typename BO, typename AggregatorFunction>
    bool operator()(const BO &bo, const AggregatorFunction &)
    {
        // Limbo calls this periodically; return true to stop optimization.
        return stop_BO_ptr->load(std::memory_order_relaxed);
    }
};

struct ChemicalSynapseParams
{
    // One-direction synapse parameters.
    // - e_syn: reversal potential
    // - g_fast/s_fast/v_fast: fast (instantaneous) sigmoid component
    // - g_slow/k1/k2/s_slow/v_slow: slow gating ODE + activation parameters
    double e_syn;
    double g_fast;
    double s_fast;
    double v_fast;
    double g_slow;
    double k1;
    double k2;
    double s_slow;
    double v_slow;
};

struct Candidate
{
    // Candidate parameters for both directions (1->2 and 2->1).
    ChemicalSynapseParams params_12;
    ChemicalSynapseParams params_21;
};

struct ChemicalSynapseEvaluation
{
    // Two-objective evaluation returned to BO.
    // - i_range_score: how close observed min/max currents are to target
    // - i_shape_score: how well current waveforms match the reference shape
    double i_range_score;
    double i_shape_score;
};

inline void copy_selected_synapse_params(ChemicalSynapseParams &runtime_params,
                                         const ChemicalSynapseParams &params,
                                         unsigned int use_i_fast,
                                         unsigned int use_i_slow)
{
    // Copy only the enabled components (fast/slow) so disabled parts remain untouched.
    if (use_i_fast || use_i_slow)
    {
        runtime_params.e_syn = params.e_syn;

        if (use_i_fast)
        {
            runtime_params.g_fast = params.g_fast;
            runtime_params.s_fast = params.s_fast;
            runtime_params.v_fast = params.v_fast;
        }

        if (use_i_slow)
        {
            runtime_params.g_slow = params.g_slow;
            runtime_params.v_slow = params.v_slow;
            runtime_params.k1 = params.k1;
            runtime_params.k2 = params.k2;
            runtime_params.s_slow = params.s_slow;
        }
    }
}

#endif
