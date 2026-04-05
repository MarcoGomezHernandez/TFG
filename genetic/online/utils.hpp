#ifndef FITNESS_H
#define FITNESS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <kfr/all.hpp>
using namespace kfr;

namespace GeneticPublicConfig
{
    inline constexpr double ETA = 0.2;
    inline constexpr double ACTIVE_WAIT_MS = 10.0;

    inline constexpr double V_MARGIN_FACTOR = 0.125;
    inline constexpr double V_FAST_MIN_FACTOR = 0.25;

    inline constexpr double EXPECTED_I_MARGIN_FACTOR = 0.5;

    inline constexpr double S_FAST_MIN_FACTOR = 4.12;
    inline constexpr double S_FAST_MAX_FACTOR = 13.72;
    inline constexpr double S_SLOW_MIN_FACTOR = 1.72;
    inline constexpr double S_SLOW_MAX_FACTOR = 3.43;

    inline constexpr double K1_MIN = 0.0000001;
    inline constexpr double K1_MAX = 1.0;
    inline constexpr double K2_MIN = 0.0000001;
    inline constexpr double K2_MAX = 1.0;

    inline constexpr double E_SYN_TERM = 3.86;

    inline constexpr double G_MIN_FACTOR = 0.0001;
}

static inline double chemical_sigmoid(double s,
                                      double v_threshold,
                                      double v_pre)
{
    return 1.0 / (1.0 + std::exp(s * (v_threshold - v_pre)));
}

struct GeneticRanges
{
    struct ParamRange
    {
        double min;
        double max;
        double mut_factor;

        ParamRange() = default;

        ParamRange(double min_value,
                   double max_value,
                   double eta = GeneticPublicConfig::ETA)
            : min(min_value),
              max(max_value),
              mut_factor(eta * (max_value - min_value))
        {
        }
    };

    ParamRange s_fast;
    ParamRange s_slow;

    double e_syn;

    ParamRange log_k1;
    ParamRange log_k2;
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
              unsigned int search_phase)
    {
        constexpr double G_MIN_FACTOR = GeneticPublicConfig::G_MIN_FACTOR;
        constexpr double K1_MAX = GeneticPublicConfig::K1_MAX;
        constexpr double K2_MIN = GeneticPublicConfig::K2_MIN;

        const double v_pre_range = v_pre_max - v_pre_min;

        const double e_syn_final_term = GeneticPublicConfig::E_SYN_TERM * (v_post_max - v_post_min);
        e_syn = search_phase ? (v_post_max + e_syn_final_term) : (v_post_min - e_syn_final_term);

        const double v_pre_margin = v_pre_range * GeneticPublicConfig::V_MARGIN_FACTOR;
        const double v_pre_margin_max = v_pre_max + v_pre_margin;
        const double v_pre_margin_min = v_pre_min - v_pre_margin;

        const double expected_i_margin = (expected_i_max - expected_i_min) * GeneticPublicConfig::EXPECTED_I_MARGIN_FACTOR;
        const double expected_i_margin_max = expected_i_max + expected_i_margin;
        const double expected_i_margin_min = expected_i_min - expected_i_margin;

        if (use_i_fast)
        {
            v_fast = ParamRange(v_pre_min + (v_pre_range * GeneticPublicConfig::V_FAST_MIN_FACTOR),
                                v_pre_margin_max);
            s_fast = ParamRange(GeneticPublicConfig::S_FAST_MIN_FACTOR / v_pre_range,
                                GeneticPublicConfig::S_FAST_MAX_FACTOR / v_pre_range);

            const double sigmoid_fast_max = chemical_sigmoid(s_fast.max, v_fast.min, v_pre_max);
            const double g_fast_max = std::max(std::abs(expected_i_margin_max / ((v_post_max - e_syn) * sigmoid_fast_max)),
                                               std::abs(expected_i_margin_min / ((v_post_min - e_syn) * sigmoid_fast_max)));
            const double g_fast_min = g_fast_max * G_MIN_FACTOR;
            log_g_fast = ParamRange(std::log(g_fast_min),
                                    std::log(g_fast_max));
        }

        if (use_i_slow)
        {
            v_slow = ParamRange(v_pre_margin_min,
                                v_pre_margin_max);
            s_slow = ParamRange(GeneticPublicConfig::S_SLOW_MIN_FACTOR / v_pre_range,
                                GeneticPublicConfig::S_SLOW_MAX_FACTOR / v_pre_range);

            log_k1 = ParamRange(std::log(GeneticPublicConfig::K1_MIN),
                                std::log(K1_MAX));
            log_k2 = ParamRange(std::log(K2_MIN),
                                std::log(GeneticPublicConfig::K2_MAX));

            const double m_max_term = K1_MAX * chemical_sigmoid(s_slow.max, v_slow.min, v_pre_max);
            const double m_max = m_max_term / (m_max_term + K2_MIN);

            const double g_slow_max = std::max(std::abs(expected_i_margin_max / ((v_post_max - e_syn) * m_max)),
                                               std::abs(expected_i_margin_min / ((v_post_min - e_syn) * m_max)));
            const double g_slow_min = g_slow_max * G_MIN_FACTOR;

            log_g_slow = ParamRange(std::log(g_slow_min),
                                    std::log(g_slow_max));
        }
    }

    GeneticRanges() = default;
};

namespace FitnessPublicConfig
{
    inline constexpr double PAD_LEN_FACTOR = 1.5;
    inline constexpr double FILTER_FC = 0.3; // En KHz, fijo, no se debe cambiar
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

struct ChemicalSynapseVariationParams
{
    double g_fast;
    double g_slow;
    double s_fast;
    double v_fast;
    double v_slow;
    double k1;
    double k2;
    double s_slow;
};

static inline void apply_variation_params(ChemicalSynapseParams &params,
                                          const ChemicalSynapseVariationParams &variation_params,
                                          unsigned int use_i_fast,
                                          unsigned int use_i_slow)
{
    if (use_i_fast)
    {
        params.g_fast = variation_params.g_fast;
        params.s_fast = variation_params.s_fast;
        params.v_fast = variation_params.v_fast;
    }

    if (use_i_slow)
    {
        params.g_slow = variation_params.g_slow;
        params.v_slow = variation_params.v_slow;
        params.k1 = variation_params.k1;
        params.k2 = variation_params.k2;
        params.s_slow = variation_params.s_slow;
    }
}

struct Individual
{
    ChemicalSynapseVariationParams variation_params_12;
    ChemicalSynapseVariationParams variation_params_21;
    double fitness;
};

#endif
