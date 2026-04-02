#ifndef FITNESS_H
#define FITNESS_H

#include <array>
#include <cmath>
#include <kfr/all.hpp>
using namespace kfr;

namespace GeneticPublicConfig
{
    inline constexpr double ETA = 0.2;

    inline constexpr double V_TH_TERM_1 = 8.0;
    inline constexpr double V_TH_TERM_2 = 4.0;

    inline constexpr double S_FAST_MIN_TERM = 4.12;
    inline constexpr double S_FAST_MAX_TERM = 13.72;
    inline constexpr double S_SLOW_MIN_TERM = 1.72;
    inline constexpr double S_SLOW_MAX_TERM = 3.43;

    inline constexpr double E_SYN_TERM = 3.86;

    inline constexpr double LOG_K_MIN = std::log(0.0000001);
    inline constexpr double LOG_K_MAX = std::log(1.0);

    inline constexpr double LOG_G_MIN = std::log(0.001);
    inline constexpr double LOG_G_MAX = std::log(5.0);
}

struct GeneticRanges
{
    double e_syn;

    double s_fast_min, s_fast_max;
    double s_slow_min, s_slow_max;

    double log_k1_min, log_k1_max;
    double log_k2_min, log_k2_max;
    double log_g_fast_min, log_g_fast_max;
    double log_g_slow_min, log_g_slow_max;

    double v_fast_min, v_fast_max;
    double v_slow_min, v_slow_max;

    double log_k1_mut_factor;
    double log_k2_mut_factor;
    double log_g_fast_mut_factor;
    double log_g_slow_mut_factor;
    double s_fast_mut_factor;
    double s_slow_mut_factor;
    double v_fast_mut_factor;
    double v_slow_mut_factor;

    init(double v_pre_min,
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
        constexpr double ETA = GeneticPublicConfig::ETA;

        const double e_syn_final_term = GeneticPublicConfig::E_SYN_TERM * (v_post_max - v_post_min);
        e_syn = search_phase ? (v_post_max + e_syn_final_term) : (v_post_min - e_syn_final_term);

        const double v_th_margin = (v_pre_range / GeneticPublicConfig::V_TH_TERM_1);
        const double v_th_max = v_pre_max + v_th_margin;
        const double v_th_min = v_pre_min - v_th_margin;

        if (use_i_fast)
        {
            v_fast_max = v_th_max;
            v_fast_min = v_th_min + ((v_th_max - v_th_min) / GeneticPublicConfig::V_TH_TERM_2);
            v_fast_mut_factor = ETA * (v_fast_max - v_fast_min);

            s_fast_max = GeneticPublicConfig::S_FAST_MAX_TERM / v_pre_range;
            s_fast_min = GeneticPublicConfig::S_FAST_MIN_TERM / v_pre_range;
            s_fast_mut_factor = ETA * (s_fast_max - s_fast_min);

            log_g_fast_min = LOG_G_MIN;
            log_g_fast_max = LOG_G_MAX;
            log_g_fast_mut_factor = ETA * (log_g_fast_max - log_g_fast_min);
        }

        if (use_i_slow)
        {
            v_slow_max = v_th_max;
            v_slow_min = v_th_min;
            v_slow_mut_factor = ETA * (v_slow_max - v_slow_min);

            s_slow_max = GeneticPublicConfig::S_SLOW_MAX_TERM / v_pre_range;
            s_slow_min = GeneticPublicConfig::S_SLOW_MIN_TERM / v_pre_range;
            s_slow_mut_factor = ETA * (s_slow_max - s_slow_min);

            log_g_slow_min = LOG_G_MIN;
            log_g_slow_max = LOG_G_MAX;
            log_g_slow_mut_factor = ETA * (log_g_slow_max - log_g_slow_min);

            log_k1_min = LOG_K_MIN;
            log_k1_max = LOG_K_MAX;
            log_k1_mut_factor = ETA * (log_k1_max - log_k1_min);
            log_k2_min = LOG_K_MIN;
            log_k2_max = LOG_K_MAX;
            log_k2_mut_factor = ETA * (log_k2_max - log_k2_min);
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
