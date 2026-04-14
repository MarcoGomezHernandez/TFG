#ifndef EVALUATION_UTILS_H
#define EVALUATION_UTILS_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <Eigen/Core>
#include <kfr/all.hpp>

namespace BOPublicConfig
{
    inline constexpr double EXPECTED_I_MARGIN_RANGE_FACTOR = 0.8;

    inline constexpr double ACTIVE_WAIT_MS = 10.0;

    inline constexpr double V_MARGIN_FACTOR = 0.125;
    inline constexpr double V_FAST_MIN_FACTOR = 0.25;

    inline constexpr double EXPECTED_I_MARGIN_FACTOR = 0.5;

    inline constexpr double S_FAST_MIN_FACTOR = 4.12;
    inline constexpr double S_FAST_MAX_FACTOR = 13.72;
    inline constexpr double S_SLOW_MIN_FACTOR = 1.72;
    inline constexpr double S_SLOW_MAX_FACTOR = 3.43;

    inline constexpr double K1_FACTOR = 3.33;
    inline constexpr double K1_MIN_FACTOR = 0.000001;

    inline constexpr double E_SYN_FAR_TERM = 3.86;
    inline constexpr double E_SYN_NEAR_TERM = 0.2;

    inline constexpr double G_MIN_FACTOR = 0.001;

    inline constexpr double R_MAX = 100.0;
    inline constexpr double R_MIN = 0.01;
}

namespace BOPublicConstants
{
    inline constexpr double SMALL_LOG = std::numeric_limits<double>::min();
    inline constexpr double SMALL_DIVISOR = std::numeric_limits<double>::epsilon();
    inline constexpr double NEGATIVE_SMALL_DIVISOR = -SMALL_DIVISOR;
}

struct StopEvaluation
{
};

inline double chemical_sigmoid(double s,
                               double v_threshold,
                               double v_pre)
{
    return 1.0 / (1.0 + std::exp(s * (v_threshold - v_pre)));
}

inline double calculate_expected_i_max_dist(double expected_i_min,
                                     double expected_i_max)
{
    const double range = expected_i_max - expected_i_min;
    return (range + (range * BOPublicConfig::EXPECTED_I_MARGIN_FACTOR * 2.0)) *
           BOPublicConfig::EXPECTED_I_MARGIN_RANGE_FACTOR;
}

inline double safe_divisor(double divisor)
{
    return std::abs(divisor) < BOPublicConstants::SMALL_DIVISOR
               ? (divisor < 0.0 ? BOPublicConstants::NEGATIVE_SMALL_DIVISOR : BOPublicConstants::SMALL_DIVISOR)
               : divisor;
}

struct BOParamRanges
{
    struct ParamRange
    {
        double min;
        double range;

        ParamRange() = default;

        ParamRange(double min, double max)
            : min(min), range(max - min)
        {
        }
    };

    ParamRange s_fast;
    ParamRange s_slow;
    ParamRange e_syn;
    ParamRange log_k1;
    ParamRange log_R;
    ParamRange log_g_fast;
    ParamRange log_g_slow;
    ParamRange v_fast;
    ParamRange v_slow;

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

        const double v_pre_range = v_pre_max - v_pre_min;
        const double v_post_range = v_post_max - v_post_min;

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

        const double v_pre_margin = v_pre_range * BOPublicConfig::V_MARGIN_FACTOR;
        const double v_pre_margin_min = v_pre_min - v_pre_margin;
        const double v_pre_margin_max = v_pre_max + v_pre_margin;

        const double expected_i_margin = (expected_i_max - expected_i_min) * BOPublicConfig::EXPECTED_I_MARGIN_FACTOR;
        const double expected_i_margin_min = expected_i_min - expected_i_margin;
        const double expected_i_margin_max = expected_i_max + expected_i_margin;
        const double safe_v_pre_range = safe_divisor(v_pre_range);

        if (use_i_fast)
        {
            v_fast = ParamRange(v_pre_min + (v_pre_range * BOPublicConfig::V_FAST_MIN_FACTOR),
                                v_pre_margin_max);
            const double s_fast_max = BOPublicConfig::S_FAST_MAX_FACTOR / safe_v_pre_range;
            s_fast = ParamRange(BOPublicConfig::S_FAST_MIN_FACTOR / safe_v_pre_range,
                                s_fast_max);

            const double sigmoid_fast_max = chemical_sigmoid(s_fast_max, v_fast.min, v_pre_max);
            const double g_fast_max = std::max(std::abs(expected_i_margin_max / safe_divisor((v_post_max - e_syn_min) * sigmoid_fast_max)),
                                               std::abs(expected_i_margin_min / safe_divisor((v_post_min - e_syn_max) * sigmoid_fast_max)));
            const double g_fast_min = g_fast_max * G_MIN_FACTOR;
            log_g_fast = ParamRange(std::log(g_fast_min == 0.0 ? SMALL_LOG : g_fast_min),
                                    std::log(g_fast_max == 0.0 ? SMALL_LOG : g_fast_max));
        }

        if (use_i_slow)
        {
            const double k1_max = BOPublicConfig::K1_FACTOR * fc;
            const double k1_min = k1_max * BOPublicConfig::K1_MIN_FACTOR;

            v_slow = ParamRange(v_pre_margin_min,
                                v_pre_margin_max);
            const double s_slow_max = BOPublicConfig::S_SLOW_MAX_FACTOR / safe_v_pre_range;
            s_slow = ParamRange(BOPublicConfig::S_SLOW_MIN_FACTOR / safe_v_pre_range,
                                s_slow_max);

            log_k1 = ParamRange(std::log(k1_min == 0.0 ? SMALL_LOG : k1_min),
                                std::log(k1_max == 0.0 ? SMALL_LOG : k1_max));

            log_R = ParamRange(std::log(R_MIN == 0.0 ? SMALL_LOG : R_MIN),
                               std::log(R_MAX == 0.0 ? SMALL_LOG : R_MAX));

            const double k2_min = k1_min * R_MIN;

            const double m_max_term = k1_max * chemical_sigmoid(s_slow_max, v_slow.min, v_pre_max);
            const double m_max = m_max_term / safe_divisor(m_max_term + k2_min);

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
    static inline std::atomic<bool> *stop_BO_ptr = nullptr;

    template <typename BO, typename AggregatorFunction>
    bool operator()(const BO &bo, const AggregatorFunction &)
    {
        return stop_BO_ptr->load(std::memory_order_relaxed);
    }
};

struct WeightedSumAggregator
{
    double w0;
    double w1;
    double total_w;

    WeightedSumAggregator(double w0 = 0.5, double w1 = 0.5)
        : w0(w0), w1(w1), total_w(w0 + w1)
    {
    }

    double operator()(const Eigen::VectorXd &y) const
    {
        assert(y.size() == 2);
        return ((w0 * y(0)) + (w1 * y(1))) / total_w;
    }
};

struct ChemicalSynapseParams
{
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
    ChemicalSynapseParams params_12;
    ChemicalSynapseParams params_21;
};

struct ChemicalSynapseEvaluation
{
    double i_range_score;
    double i_shape_score;
};

inline void copy_selected_synapse_params(ChemicalSynapseParams &runtime_params,
                                         const ChemicalSynapseParams &params,
                                         unsigned int use_i_fast,
                                         unsigned int use_i_slow)
{
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
