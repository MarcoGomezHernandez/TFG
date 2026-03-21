
#include <algorithm>
#include <cmath>
#include <iostream>
#include <kfr/all.hpp>
using namespace kfr;

namespace FitnessConfig
{
    static constexpr double SYN_RANGE_WEIGHT = 0.5;
    static constexpr double VPRE_SYN_COMP_WEIGHT = 0.5;

    static constexpr double IFAST_WEIGHT = 0.5;
    static constexpr double ISLOW_WEIGHT = 0.5;

    static constexpr size_t FILTER_PAD_LEN = 1000;
    static constexpr int BUTTERWORTH_ORDER = 4;
}

static inline double rescale_to_target(double value,
                                       double src_min, double src_max,
                                       double dst_min, double dst_max)
{
    const double norm = (value - src_min) / (src_max - src_min);
    return dst_min + (norm * (dst_max - dst_min));
}

static inline double pearson_score(univector<double> &sig,
                                   const univector<double> &ref_sig_centered,
                                   double ref_sig_factor,
                                   bool search_phase)
{
    sig -= mean(sig);
    const double sig_factor = std::sqrt(sum(sqr(sig)));
    const double r = sum(sig * ref_sig_centered) / (sig_factor * ref_sig_factor);
    const double normalized = (r + 1.0) / 2.0;
    return search_phase ? normalized : 1.0 - normalized;
}

static inline double calc_range_score_component(double observed_min, double observed_max,
                                                double expected_min, double expected_max,
                                                double max_diff)
{
    return 1.0 - (((std::abs(observed_min - expected_min) +
                    std::abs(observed_max - expected_max)) *
                   0.5) /
                  max_diff);
}


static inline double score_one_direction(
    const univector_ref<const double> &vpre_sig,
    const univector_ref<const double> &ref_ifast_sig,
    const univector_ref<const double> &ref_islow_sig,
    univector<double> &cand_ifast_sig,
    univector<double> &cand_islow_sig,
    double expected_i_min,
    double expected_i_max,
    bool search_phase,
    bool use_ifast,
    bool use_islow,
    univector<double> &ref_ifast_centered,
    univector<double> &ref_islow_centered)
{
    constexpr double IFAST_WEIGHT = FitnessConfig::IFAST_WEIGHT;
    constexpr double ISLOW_WEIGHT = FitnessConfig::ISLOW_WEIGHT;

    const bool use_both = use_ifast && use_islow;

    double vpre_syn_comp_accum = 0.0;
    double i_range_accum = 0.0;
    double total_weight = 0.0;

    if (use_ifast)
    {
        ref_ifast_centered = ref_ifast_sig - mean(ref_ifast_sig);
        const double ref_ifast_factor = std::sqrt(sum(sqr(ref_ifast_centered)));

        double ifast_min_to_fit, ifast_max_to_fit;
        if (use_both)
        {
            const double vpre_min = minof(vpre_sig);
            const double vpre_max = maxof(vpre_sig);
            ifast_min_to_fit = rescale_to_target(minof(ref_ifast_sig), vpre_min, vpre_max,
                                                 expected_i_min, expected_i_max);
            ifast_max_to_fit = rescale_to_target(maxof(ref_ifast_sig), vpre_min, vpre_max,
                                                 expected_i_min, expected_i_max);
        }
        else
        {
            ifast_min_to_fit = expected_i_min;
            ifast_max_to_fit = expected_i_max;
        }
        const double ifast_range_to_fit = ifast_max_to_fit - ifast_min_to_fit;

        const double vpost_ifast_comp_score =
            pearson_score(cand_ifast_sig,
                          ref_ifast_centered,
                          ref_ifast_factor,
                          search_phase);
        vpre_syn_comp_accum += IFAST_WEIGHT * vpost_ifast_comp_score;

        const double comp_ifast = calc_range_score_component(
            minof(cand_ifast_sig),
            maxof(cand_ifast_sig),
            ifast_min_to_fit,
            ifast_max_to_fit,
            ifast_range_to_fit);
        i_range_accum += IFAST_WEIGHT * comp_ifast;

        total_weight += IFAST_WEIGHT;
    }

    if (use_islow)
    {
        ref_islow_centered = ref_islow_sig - mean(ref_islow_sig);
        const double ref_islow_factor = std::sqrt(sum(sqr(ref_islow_centered)));

        double islow_min_to_fit, islow_max_to_fit;
        if (use_both)
        {
            const double vpre_min = minof(vpre_sig);
            const double vpre_max = maxof(vpre_sig);
            islow_min_to_fit = rescale_to_target(minof(ref_islow_sig), vpre_min, vpre_max,
                                                 expected_i_min, expected_i_max);
            islow_max_to_fit = rescale_to_target(maxof(ref_islow_sig), vpre_min, vpre_max,
                                                 expected_i_min, expected_i_max);
        }
        else
        {
            islow_min_to_fit = expected_i_min;
            islow_max_to_fit = expected_i_max;
        }
        const double islow_range_to_fit = islow_max_to_fit - islow_min_to_fit;

        const double vpost_islow_comp_score =
            pearson_score(cand_islow_sig,
                          ref_islow_centered,
                          ref_islow_factor,
                          search_phase);
        vpre_syn_comp_accum += ISLOW_WEIGHT * vpost_islow_comp_score;

        const double comp_islow = calc_range_score_component(
            minof(cand_islow_sig),
            maxof(cand_islow_sig),
            islow_min_to_fit,
            islow_max_to_fit,
            islow_range_to_fit);
        i_range_accum += ISLOW_WEIGHT * comp_islow;

        total_weight += ISLOW_WEIGHT;
    }

    const double final_score = ((FitnessConfig::SYN_RANGE_WEIGHT * i_range_accum) +
                                (FitnessConfig::VPRE_SYN_COMP_WEIGHT * vpre_syn_comp_accum)) /
                               total_weight;

    return final_score;
}

double fitness_from_sigs(const FitnessExtraData &extra,
                         SigBuffers &buffers)
{

    univector<double> &v1_sig = buffers.v1_sig;
    univector<double> &v2_sig = buffers.v2_sig;
    univector<double> &syn1_ifast_sig = buffers.syn1_ifast_sig;
    univector<double> &syn1_islow_sig = buffers.syn1_islow_sig;
    univector<double> &syn2_ifast_sig = buffers.syn2_ifast_sig;
    univector<double> &syn2_islow_sig = buffers.syn2_islow_sig;

    univector<double> &ref_ifast_centered = buffers.ref_ifast_centered;
    univector<double> &ref_islow_centered = buffers.ref_islow_centered;

    const bool use_ifast = extra.use_ifast;
    const bool use_islow = extra.use_islow;
    const bool search_phase = extra.phase;


    const double score_dir1 = score_one_direction(
        v1_sig,
        extra.ifast_sig_1, extra.islow_sig_1,
        syn1_ifast_sig, syn1_islow_sig,
        extra.i_range_expected_min_syn1, extra.i_range_expected_max_syn1,
        search_phase,
        use_ifast, use_islow,
        ref_ifast_centered, ref_islow_centered);


    const double score_dir2 = score_one_direction(
        v2_sig,
        extra.ifast_sig_2, extra.islow_sig_2,
        syn2_ifast_sig, syn2_islow_sig,
        extra.i_range_expected_min_syn2, extra.i_range_expected_max_syn2,
        search_phase,
        use_ifast, use_islow,
        ref_ifast_centered, ref_islow_centered);


    const double final_score = (score_dir1 + score_dir2) / 2.0;

    if (!(std::isfinite(final_score)) || final_score < 0.0)
        return 0.0;

    return final_score;
}
