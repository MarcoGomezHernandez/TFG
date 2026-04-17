#ifndef EVALUATION_H
#define EVALUATION_H

#include <cmath>
#include <cstddef>
#include <kfr/all.hpp>
#include <ChemicalSynapsis.h>
#include "scaling.hpp"
#include "utils.hpp"

// Public constants for candidate-scoring helpers.
namespace EvaluationPublicConfig
{
    inline constexpr double VERY_BAD_RANGE_SCORE = -1e6;
}

struct ChemicalSynapseEvaluation
{
    // Score for matching expected current min/max bounds.
    double i_range_score;
    // Score for matching current waveform shape.
    double i_shape_score;

    ChemicalSynapseEvaluation(double i_range_score_ = 0.0,
                              double i_shape_score_ = 0.0)
        : i_range_score(std::isfinite(i_range_score_) ? i_range_score_ : EvaluationPublicConfig::VERY_BAD_RANGE_SCORE),
          i_shape_score(std::isfinite(i_shape_score_) ? i_shape_score_ : 0.0)
    {
    }
};

struct ConstantEvaluationVals
{
    // Centered reference signals + normalization factors for Pearson score.
    kfr::univector<double> ref_i_fast_sig_centered;
    double ref_i_fast_sig_factor;
    kfr::univector<double> ref_i_slow_sig_centered;
    double ref_i_slow_sig_factor;

    // Expected ranges and normalization distances for range score.
    double ref_i_fast_min;
    double ref_i_fast_max;
    double i_fast_dist_max;

    double ref_i_slow_min;
    double ref_i_slow_max;
    double i_slow_dist_max;
};

struct EvaluationISigBuffers
{
    // Allocate only enabled component buffers to avoid unnecessary memory use.
    EvaluationISigBuffers(size_t size_to_reserve,
                          bool use_i_fast,
                          bool use_i_slow)
    {
        if (use_i_fast)
            i_fast_sig.resize(size_to_reserve);
        if (use_i_slow)
            i_slow_sig.resize(size_to_reserve);
    }

    kfr::univector<double> i_fast_sig;
    kfr::univector<double> i_slow_sig;
};

// Precompute reference constants used by all candidate evaluations.
ConstantEvaluationVals calc_constant_evaluation_vals(
    const kfr::univector_ref<double> &v_pre_sig,
    double v_pre_min,
    double v_pre_max,
    double csv_step,
    double fc,
    bool use_i_fast,
    bool use_i_slow,
    double expected_i_min,
    double expected_i_max,
    bool search_phase);

template <typename Integrator, typename NeuronType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
// Simulate a candidate in closed-loop and compute (range_score, shape_score).
ChemicalSynapseEvaluation evaluate_candidate(
    const ChemicalSynapseParams &candidate,
    ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapse,
    NeuronType &model_neur,
    const ScaledSigResult &scaled_result,
    EvaluationISigBuffers &buffers,
    size_t v_pre_sig_start_idx,
    const ConstantEvaluationVals &constant_evaluation_vals,
    bool use_i_fast,
    bool use_i_slow,
    bool search_phase,
    double i_min,
    double i_max,
    ResetStateFuncType reset_state_neur,
    GetVFuncType get_v_neur);

#include "evaluation.tpp"

#endif
