#ifndef FITNESS_H
#define FITNESS_H

#include <array>
#include <cmath>
#include <kfr/all.hpp>
using namespace kfr;

namespace GeneticPublicConfig
{
    inline constexpr double ETA = 0.2;
    inline constexpr double ACTIVE_WAIT_SECS = 0.01;

    inline constexpr double V_TH_TERM_1 = 8.0;
    inline constexpr double V_TH_TERM_2 = 4.0;

    inline constexpr double S_FAST_MIN_TERM = 4.12;
    inline constexpr double S_FAST_MAX_TERM = 13.72;
    inline constexpr double S_SLOW_MIN_TERM = 1.72;
    inline constexpr double S_SLOW_MAX_TERM = 3.43;

    inline constexpr double LOG_K_MIN = std::log(0.0000001);
    inline constexpr double LOG_K_MAX = std::log(1.0);

    inline constexpr double LOG_G_MIN = std::log(0.001);
    inline constexpr double LOG_G_MAX = std::log(5.0);
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
              double dt,
              unsigned int use_i_fast,
              unsigned int use_i_slow,
              unsigned int search_phase)
    {
        const double v_pre_range = v_pre_max - v_pre_min;

        constexpr double LOG_G_MIN = GeneticPublicConfig::LOG_G_MIN;
        constexpr double LOG_G_MAX = GeneticPublicConfig::LOG_G_MAX;
        constexpr double LOG_K_MIN = GeneticPublicConfig::LOG_K_MIN;
        constexpr double LOG_K_MAX = GeneticPublicConfig::LOG_K_MAX;

        const double v_th_margin = (v_pre_range / GeneticPublicConfig::V_TH_TERM_1);
        const double v_th_max = v_pre_max + v_th_margin;
        const double v_th_min = v_pre_min - v_th_margin;

        if (use_i_fast)
        {
            v_fast = ParamRange(v_th_min + ((v_th_max - v_th_min) / GeneticPublicConfig::V_TH_TERM_2),
                                v_th_max);
            s_fast = ParamRange(GeneticPublicConfig::S_FAST_MIN_TERM / v_pre_range,
                                GeneticPublicConfig::S_FAST_MAX_TERM / v_pre_range);
            log_g_fast = ParamRange(LOG_G_MIN,
                                    LOG_G_MAX);
        }

        if (use_i_slow)
        {
            v_slow = ParamRange(v_th_min,
                                v_th_max);
            s_slow = ParamRange(GeneticPublicConfig::S_SLOW_MIN_TERM / v_pre_range,
                                GeneticPublicConfig::S_SLOW_MAX_TERM / v_pre_range);
            log_g_slow = ParamRange(LOG_G_MIN,
                                    LOG_G_MAX);
            log_k1 = ParamRange(LOG_K_MIN,
                                LOG_K_MAX);
            log_k2 = ParamRange(LOG_K_MIN,
                                LOG_K_MAX);
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
