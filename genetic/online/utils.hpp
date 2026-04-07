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

    inline constexpr double K1_FACTOR = 3.33;
    inline constexpr double K2_FACTOR = 3.33;

    inline constexpr double E_SYN_FAR_TERM = 3.86;
    inline constexpr double E_SYN_NEAR_TERM = 0.2;

    inline constexpr double G_MIN_FACTOR = 0.0001;
    inline constexpr double K1_MIN_FACTOR = 0.0000001;
    inline constexpr double K2_MIN_FACTOR = 0.0000001;
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

    ParamRange e_syn;

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
              unsigned int search_phase,
              double fc)
    {
        constexpr double G_MIN_FACTOR = GeneticPublicConfig::G_MIN_FACTOR;

        const double v_pre_range = v_pre_max - v_pre_min;
        const double v_post_range = v_post_max - v_post_min;

        const double e_syn_far_final_term = v_post_range * GeneticPublicConfig::E_SYN_FAR_TERM;
        const double e_syn_near_final_term = v_post_range * GeneticPublicConfig::E_SYN_NEAR_TERM;
        double e_syn_max = 0.0, e_syn_min = 0.0;
        if (search_phase)
        {
            e_syn_max = v_post_max + e_syn_far_final_term;
            e_syn_min = v_post_max + e_syn_near_final_term;
        }
        else
        {
            e_syn_max = v_post_min - e_syn_near_final_term;
            e_syn_min = v_post_min - e_syn_far_final_term;
        }
        e_syn = ParamRange(e_syn_min, e_syn_max);

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
            const double g_fast_max = std::max(std::abs(expected_i_margin_max / ((v_post_max - e_syn_min) * sigmoid_fast_max)),
                                               std::abs(expected_i_margin_min / ((v_post_min - e_syn_max) * sigmoid_fast_max)));
            const double g_fast_min = g_fast_max * G_MIN_FACTOR;
            log_g_fast = ParamRange(std::log(g_fast_min),
                                    std::log(g_fast_max));
        }

        if (use_i_slow)
        {
            const double k1_max = GeneticPublicConfig::K1_FACTOR * fc;
            const double k1_min = k1_max * GeneticPublicConfig::K1_MIN_FACTOR;
            const double k2_max = GeneticPublicConfig::K2_FACTOR * fc;
            const double k2_min = k2_max * GeneticPublicConfig::K2_MIN_FACTOR;

            v_slow = ParamRange(v_pre_margin_min,
                                v_pre_margin_max);
            s_slow = ParamRange(GeneticPublicConfig::S_SLOW_MIN_FACTOR / v_pre_range,
                                GeneticPublicConfig::S_SLOW_MAX_FACTOR / v_pre_range);

            log_k1 = ParamRange(std::log(k1_min),
                                std::log(k1_max));
            log_k2 = ParamRange(std::log(k2_min),
                                std::log(k2_max));

            const double m_max_term = k1_max * chemical_sigmoid(s_slow.max, v_slow.min, v_pre_max);
            const double m_max = m_max_term / (m_max_term + k2_min);

            const double g_slow_max = std::max(std::abs(expected_i_margin_max / ((v_post_max - e_syn_min) * m_max)),
                                               std::abs(expected_i_margin_min / ((v_post_min - e_syn_max) * m_max)));
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
    inline constexpr double FILTER_FC = 0.3; // En KHz, valor por defecto, ajustado para neurona viva
};

struct FitnessPadBuffers
{
    univector<double> padded_buff_12;
    univector<double> padded_buff_21;

    FitnessPadBuffers(size_t padded_buff_size_12,
                      size_t padded_buff_size_21)
        : padded_buff_12(padded_buff_size_12),
          padded_buff_21(padded_buff_size_21)
    {
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

struct IndividualChemicalSynapseParams
{
    double e_syn;
    double log_g_fast;
    double s_fast;
    double v_fast;
    double log_g_slow;
    double log_k1;
    double log_k2;
    double s_slow;
    double v_slow;
};

static inline void copy_individual_synapse_params_to_runtime(ChemicalSynapseParams &runtime_params,
                                                             const IndividualChemicalSynapseParams &individual_params,
                                                             unsigned int use_i_fast,
                                                             unsigned int use_i_slow)
{
    if (use_i_fast || use_i_slow)
    {
        runtime_params.e_syn = individual_params.e_syn;
        if (use_i_fast)
        {
            runtime_params.g_fast = std::exp(individual_params.log_g_fast);
            runtime_params.s_fast = individual_params.s_fast;
            runtime_params.v_fast = individual_params.v_fast;
        }

        if (use_i_slow)
        {
            runtime_params.g_slow = std::exp(individual_params.log_g_slow);
            runtime_params.v_slow = individual_params.v_slow;
            runtime_params.k1 = std::exp(individual_params.log_k1);
            runtime_params.k2 = std::exp(individual_params.log_k2);
            runtime_params.s_slow = individual_params.s_slow;
        }
    }
}

struct Individual
{
    IndividualChemicalSynapseParams params_12;
    IndividualChemicalSynapseParams params_21;
    double fitness;
};

#endif
