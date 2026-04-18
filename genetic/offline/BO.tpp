#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <kfr/all.hpp>
#include "scaling.hpp"
#include <ChemicalSynapsis.h>

using namespace limbo;

// Offline BO implementation for one synapse direction.
// It optimizes a normalized vector x in [0,1]^d and decodes x into
// physical synapse parameters before evaluating each candidate by simulation.

namespace BOPrivateConfig
{
    // Limbo hyper-parameter optimization period bounds.
    static constexpr int HP_PERIOD_MIN = 1;
    static constexpr int HP_PERIOD_DIVISOR = 25;
    // Internal optimizer iteration multipliers (scaled by BO dimension).
    static constexpr int RPROP_ITER_FACTOR = 50;
    static constexpr int NLOPT_ITER_FACTOR = 150;

    // Evaluation objective weights.
    static constexpr double I_SHAPE_WEIGHT = 0.5;
    static constexpr double I_RANGE_WEIGH = 0.5;

    // Parameter-range initialization factors.
    static constexpr double V_FAST_MIN_FACTOR = 0.25;
    static constexpr double S_FAST_MIN_FACTOR = 4.12;
    static constexpr double S_FAST_MAX_FACTOR = 13.72;
    static constexpr double S_SLOW_MIN_FACTOR = 1.72;
    static constexpr double S_SLOW_MAX_FACTOR = 3.43;
    static constexpr double K1_MAX_FACTOR = 3.33;
    static constexpr double K1_MIN_FACTOR = 3.33e-6;
    static constexpr double E_SYN_FAR_TERM = 3.86;
    static constexpr double E_SYN_NEAR_TERM = 0.2;
    static constexpr double G_MIN_FACTOR = 0.001;
    static constexpr double R_MAX = 40.0;
    static constexpr double R_MIN = 0.01;
}

namespace BOPrivateConstants
{
    // Normalization denominator for weighted objective.
    static constexpr double TOTAL_WEIGHT = BOPrivateConfig::I_SHAPE_WEIGHT + BOPrivateConfig::I_RANGE_WEIGH;

    // Guard value to avoid log(0).
    static constexpr double SMALL_LOG = std::numeric_limits<double>::min();
}

struct Params
{
    struct bayes_opt_bobase : public defaults::bayes_opt_bobase
    {
        // Constrain optimizer to [0,1]^d (decoded later to physical units).
        BO_PARAM(bool, bounded, true);
        // Enable Limbo runtime statistics.
        BO_PARAM(bool, stats_enabled, true);
    };

    struct bayes_opt_boptimizer : public defaults::bayes_opt_boptimizer
    {
        // Frequency (iterations) of GP hyper-parameter updates.
        BO_DYN_PARAM(int, hp_period);
    };

    struct init_lhs : public defaults::init_lhs
    {
        // Number of Latin Hypercube initialization samples.
        BO_DYN_PARAM(int, samples);
    };

    struct kernel : public defaults::kernel
    {
        // Assume almost noiseless observations in offline deterministic evaluations.
        BO_PARAM(double, noise, 1e-6);
        BO_PARAM(bool, optimize_noise, false);
    };

    struct kernel_squared_exp_ard : public defaults::kernel_squared_exp_ard
    {
        // Number of columns in the Λ matrix for characteristic length-scales.
        // A value of 0 disables the Λ matrix, reverting to the standard ARD kernel.
        BO_PARAM(int, k, 0);
        // Initial signal variance for the SE-ARD kernel.
        BO_PARAM(double, sigma_sq, 1.0);
    };

    struct acqui_ei : public defaults::acqui_ei
    {
        // Small jitter to keep EI explotation numerically robust.
        BO_PARAM(double, jitter, 0.003);
    };

    struct opt_nloptnograd : public defaults::opt_nloptnograd
    {
        // Iteration budget for acquisition optimization.
        BO_DYN_PARAM(int, iterations);
        // Disable function change tolerance-based early stopping.
        BO_PARAM(double, fun_tolerance, -1);
        // Disable relative x-change tolerance-based early stopping.
        BO_PARAM(double, xrel_tolerance, -1);
    };

    struct opt_rprop : public defaults::opt_rprop
    {
        // Iteration budget for GP hyper-parameter optimization.
        BO_DYN_PARAM(int, iterations);
        // Termination tolerance for the Rprop optimizer.
        BO_PARAM(double, eps_stop, 1e-6);
    };

    struct stop_maxiterations : public defaults::stop_maxiterations
    {
        // Global BO iteration cap for the entire optimization loop.
        BO_DYN_PARAM(int, iterations);
    };

    struct mean_constant : public defaults::mean_constant
    {
        // Constant prior mean value for the GP.
        BO_PARAM(double, constant, 0.5);
    };
};

static double chemical_sigmoid(double s,
                               double v_threshold,
                               double v_pre)
{
    // Logistic activation used by chemical synapse equations.
    return 1.0 / (1.0 + std::exp(s * (v_threshold - v_pre)));
}

struct BOParamRanges
{
    struct ParamRange
    {
        // Decoding helper: decoded_value = min + x * range, with x in [0,1].
        double min;
        double range;

        ParamRange() = default;

        ParamRange(double min, double max)
            : min(min), range(max - min)
        {
        }
    };

    // Decoding ranges (log_* are optimized in log-domain for conditioning).
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
              bool use_i_fast,
              bool use_i_slow,
              bool search_phase,
              double fc)
    {
        // Build parameter decoding ranges from voltage/current constraints.
        constexpr double G_MIN_FACTOR = BOPrivateConfig::G_MIN_FACTOR;
        constexpr double R_MIN = BOPrivateConfig::R_MIN;
        constexpr double SMALL_LOG = BOPrivateConstants::SMALL_LOG;
        constexpr double R_MAX = BOPrivateConfig::R_MAX;

        const double v_pre_range = v_pre_max - v_pre_min;
        const double v_post_range = v_post_max - v_post_min;

        // E_syn range depends on phase search direction relative to post voltage.
        const double e_syn_far_final_term = v_post_range * BOPrivateConfig::E_SYN_FAR_TERM;
        const double e_syn_near_final_term = v_post_range * BOPrivateConfig::E_SYN_NEAR_TERM;
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

        // Expand expected current range with a margin to avoid over-constraining BO.
        const double expected_i_margin = (expected_i_max - expected_i_min) * BOPublicConfig::EXPECTED_I_MARGIN_FACTOR;
        const double expected_i_margin_min = expected_i_min - expected_i_margin;
        const double expected_i_margin_max = expected_i_max + expected_i_margin;
        const double safe_v_pre_range = safe_divisor(v_pre_range);

        if (use_i_fast)
        {
            // Fast subset bounds.
            v_fast = ParamRange(v_pre_min + (v_pre_range * BOPrivateConfig::V_FAST_MIN_FACTOR),
                                v_pre_max);
            const double s_fast_max = BOPrivateConfig::S_FAST_MAX_FACTOR / safe_v_pre_range;
            s_fast = ParamRange(BOPrivateConfig::S_FAST_MIN_FACTOR / safe_v_pre_range,
                                s_fast_max);

            const double sigmoid_fast_max = chemical_sigmoid(s_fast_max, v_fast.min, v_pre_max);
            // Bound g_fast so expected current interval remains reachable at extremes.
            const double g_fast_max = std::max(std::abs(expected_i_margin_max / safe_divisor((v_post_max - e_syn_min) * sigmoid_fast_max)),
                                               std::abs(expected_i_margin_min / safe_divisor((v_post_min - e_syn_max) * sigmoid_fast_max)));
            const double g_fast_min = g_fast_max * G_MIN_FACTOR;
            // Conductance optimized in log-domain for numerical stability.
            log_g_fast = ParamRange(std::log(g_fast_min == 0.0 ? SMALL_LOG : g_fast_min),
                                    std::log(g_fast_max == 0.0 ? SMALL_LOG : g_fast_max));
        }

        if (use_i_slow)
        {
            // Slow subset bounds.
            const double k1_max = BOPrivateConfig::K1_MAX_FACTOR * fc;
            const double k1_min = BOPrivateConfig::K1_MIN_FACTOR * fc;

            v_slow = ParamRange(v_pre_min,
                                v_pre_max);
            const double s_slow_max = BOPrivateConfig::S_SLOW_MAX_FACTOR / safe_v_pre_range;
            s_slow = ParamRange(BOPrivateConfig::S_SLOW_MIN_FACTOR / safe_v_pre_range,
                                s_slow_max);

            // k1 and R = k2/k1 are optimized in log-domain.
            log_k1 = ParamRange(std::log(k1_min == 0.0 ? SMALL_LOG : k1_min),
                                std::log(k1_max == 0.0 ? SMALL_LOG : k1_max));

            log_R = ParamRange(std::log(R_MIN == 0.0 ? SMALL_LOG : R_MIN),
                               std::log(R_MAX == 0.0 ? SMALL_LOG : R_MAX));

            const double k2_min = k1_min * R_MIN;

            const double m_max_term = k1_max * chemical_sigmoid(s_slow_max, v_slow.min, v_pre_max);
            // Conservative upper bound of steady-state m used to bound g_slow.
            const double m_max = m_max_term / safe_divisor(m_max_term + k2_min);

            // Bound g_slow so expected current interval remains reachable at extremes.
            const double g_slow_max = std::max(std::abs(expected_i_margin_max / safe_divisor((v_post_max - e_syn_min) * m_max)),
                                               std::abs(expected_i_margin_min / safe_divisor((v_post_min - e_syn_max) * m_max)));
            const double g_slow_min = g_slow_max * G_MIN_FACTOR;

            log_g_slow = ParamRange(std::log(g_slow_min == 0.0 ? SMALL_LOG : g_slow_min),
                                    std::log(g_slow_max == 0.0 ? SMALL_LOG : g_slow_max));
        }
    }

    BOParamRanges() = default;
};

using Model_t = model::GP<
    Params,
    kernel::SquaredExpARD<Params>,
    mean::Constant<Params>,
    model::gp::KernelLFOpt<Params, opt::Rprop<Params>>>;

// Expected-Improvement acquisition over the GP posterior.
using Acqui_t = acqui::EI<Params, Model_t>;

// Gradient-free acquisition optimizer (LN_SBPLX is a local, derivative-free method suitable for low-dimensional problems).
using AcquiOpt_t = opt::NLOptNoGrad<Params, nlopt::LN_SBPLX>;

using Stat_t = boost::fusion::vector<
    stat::Samples<Params>,
    stat::Observations<Params>,
    stat::GPAcquisitions<Params>,
    stat::BestObservations<Params>,
    stat::AggregatedObservations<Params>>;

using Init_t = init::LHS<Params>;

using Stop_t = stop::MaxIterations<Params>;

using BO_t = bayes_opt::BOptimizer<
    Params,
    stopcrit<Stop_t>,
    modelfun<Model_t>,
    acquifun<Acqui_t>,
    statsfun<Stat_t>,
    initfun<Init_t>,
    acquiopt<AcquiOpt_t>>;

static double decode_param(double x_val, const BOParamRanges::ParamRange &range)
{
    // Clamp to [0,1] then map to physical interval.
    return range.min + (std::clamp(x_val, 0.0, 1.0) * range.range);
}

static ChemicalSynapseParams decode_to_candidate(const Eigen::VectorXd &x,
                                                 const BOParamRanges &ranges,
                                                 bool use_i_fast,
                                                 bool use_i_slow)
{
    // Decode normalized BO vector into one-direction synapse parameters.
    ChemicalSynapseParams candidate{};
    size_t idx = 0;

    candidate.e_syn = decode_param(x(idx++), ranges.e_syn);

    if (use_i_fast)
    {
        // Conductance optimized in log-domain -> exp back to linear.
        candidate.g_fast = std::exp(decode_param(x(idx++), ranges.log_g_fast));
        candidate.s_fast = decode_param(x(idx++), ranges.s_fast);
        candidate.v_fast = decode_param(x(idx++), ranges.v_fast);
    }

    if (use_i_slow)
    {
        candidate.g_slow = std::exp(decode_param(x(idx++), ranges.log_g_slow));
        candidate.v_slow = decode_param(x(idx++), ranges.v_slow);
        candidate.k1 = std::exp(decode_param(x(idx++), ranges.log_k1));
        // R = k2/k1 is optimized to avoid strongly coupled raw k1/k2 search.
        const double R = std::exp(decode_param(x(idx++), ranges.log_R));
        candidate.k2 = candidate.k1 * R;
        candidate.s_slow = decode_param(x(idx++), ranges.s_slow);
    }

    return candidate;
}

struct EvaluationFunctorBase
{
    // Limbo dynamic parameter: BO input dimensionality.
    BO_DYN_PARAM(int, dim_in);
};

template <typename Integrator, typename NeuronType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
struct EvaluationFunctor : public EvaluationFunctorBase
{
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    // Adapter between Limbo's optimizer interface and the domain evaluator.
    // Holds references to simulation objects and precomputed scoring constants.

    EvaluationFunctor(ChemicalSynapsisType &synapse_,
                      NeuronType &model_neur_,
                      const ScaledSigResult &scaled_result_,
                      EvaluationISigBuffers &buffers_,
                      const ConstantEvaluationVals &constant_evaluation_vals_,
                      const BOParamRanges &ranges_,
                      bool use_i_fast_,
                      bool use_i_slow_,
                      bool search_phase_,
                      size_t v_pre_sig_start_idx_,
                      double i_min_,
                      double i_max_,
                      bool verbose_,
                      ResetStateFuncType reset_state_neur_,
                      GetVFuncType get_v_neur_)
        : synapse(synapse_),
          model_neur(model_neur_),
          scaled_result(scaled_result_),
          buffers(buffers_),
          constant_evaluation_vals(constant_evaluation_vals_),
          ranges(ranges_),
          use_i_fast(use_i_fast_),
          use_i_slow(use_i_slow_),
          search_phase(search_phase_),
          v_pre_sig_start_idx(v_pre_sig_start_idx_),
          i_min(i_min_),
          i_max(i_max_),
          verbose(verbose_),
          reset_state_neur(reset_state_neur_),
          get_v_neur(get_v_neur_)
    {
        // Plain aggregate constructor used by Limbo callback wrapper.
    }

    ChemicalSynapsisType &synapse;
    NeuronType &model_neur;
    const ScaledSigResult &scaled_result;
    EvaluationISigBuffers &buffers;
    const ConstantEvaluationVals &constant_evaluation_vals;
    const BOParamRanges &ranges;
    bool use_i_fast;
    bool use_i_slow;
    bool search_phase;
    size_t v_pre_sig_start_idx;
    double i_min;
    double i_max;
    bool verbose;
    ResetStateFuncType reset_state_neur;
    GetVFuncType get_v_neur;
    mutable size_t evaluation_count = 0;

    BO_PARAM(int, dim_out, 1);

    Eigen::VectorXd operator()(const Eigen::VectorXd &x) const
    {
        // Evaluate one BO sample and return scalar objective as 1D vector.
        ChemicalSynapseParams candidate = decode_to_candidate(x, ranges, use_i_fast, use_i_slow);

        ChemicalSynapseEvaluation evaluations = evaluate_candidate(
            candidate,
            synapse,
            model_neur,
            scaled_result,
            buffers,
            v_pre_sig_start_idx,
            constant_evaluation_vals,
            use_i_fast,
            use_i_slow,
            search_phase,
            i_min,
            i_max,
            reset_state_neur,
            get_v_neur);

        const double y = ((evaluations.i_range_score * BOPrivateConfig::I_RANGE_WEIGH) +
                          (evaluations.i_shape_score * BOPrivateConfig::I_SHAPE_WEIGHT)) /
                         BOPrivateConstants::TOTAL_WEIGHT;
        // Limbo maximizes this scalar objective.

        if (verbose)
        {
            // Track per-evaluation objective decomposition.
            std::cout << "Evaluation " << ++evaluation_count << ": " << y << " (range: " << evaluations.i_range_score << ", shape: " << evaluations.i_shape_score << ")" << std::endl;
        }

        return limbo::tools::make_vector(y);
    }
};

BO_DECLARE_DYN_PARAM(int, EvaluationFunctorBase, dim_in);
BO_DECLARE_DYN_PARAM(int, Params::bayes_opt_boptimizer, hp_period);
BO_DECLARE_DYN_PARAM(int, Params::init_lhs, samples);
BO_DECLARE_DYN_PARAM(int, Params::stop_maxiterations, iterations);
BO_DECLARE_DYN_PARAM(int, Params::opt_rprop, iterations);
BO_DECLARE_DYN_PARAM(int, Params::opt_nloptnograd, iterations);

static inline size_t count_params(bool use_i_fast,
                                  bool use_i_slow)
{
    // e_syn is always optimized; fast contributes 3 vars, slow contributes 5.
    size_t num_params = 1;
    if (use_i_fast)
        num_params += 3;
    if (use_i_slow)
        num_params += 5;
    return num_params;
}

template <typename Integrator, typename NeuronType,
          CreateFunc<NeuronType> CreateFuncType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
std::optional<ChemicalSynapseParams> BO(const std::string &csv_path,
                                        size_t column_idx,
                                        double csv_step,
                                        double start_time,
                                        double stabilization_time,
                                        double evaluation_time,
                                        double observation_time,
                                        size_t initial_samples,
                                        size_t iterations,
                                        NumericIntegrator integrator,
                                        NeuronModel model,
                                        bool search_phase,
                                        bool check_drift,
                                        SynComponent syn_component,
                                        CreateFuncType create_neur,
                                        ResetStateFuncType reset_state_neur,
                                        GetVFuncType get_v_neur,
                                        typename NeuronType::variable neur_v_var,
                                        int syn_model_step_factor,
                                        double fc,
                                        double expected_i_min,
                                        double expected_i_max,
                                        double i_min,
                                        double i_max,
                                        bool verbose)
{
    // Offline BO orchestration pipeline:
    // 1) Load and scale external voltage trace.
    // 2) Build parameter decoding ranges from observed signal + constraints.
    // 3) Simulate candidates and compute range/shape score.
    // 4) Return best decoded synapse parameters.
    if (csv_step <= 0.0 || evaluation_time <= 0.0 || observation_time <= 0.0 || stabilization_time < 0.0 || fc <= 0.0 || expected_i_min >= expected_i_max || i_min > i_max)
    {
        throw std::invalid_argument("Invalid arguments: csv_step, evaluation_time, observation_time, and fc must be positive; stabilization_time non-negative; expected_i_min must be less than expected_i_max; i_min must be less or equal to i_max");
    }

    const bool use_i_fast = (syn_component != SynComponent::ISLOW);
    const bool use_i_slow = (syn_component != SynComponent::IFAST);

    if (!use_i_fast && !use_i_slow)
    {
        throw std::invalid_argument("At least one of use_i_fast or use_i_slow must be true");
    }

    const std::optional<ScaledSigResult> scaled_result_opt = scale_sig(
        csv_path, column_idx, csv_step, start_time, evaluation_time + stabilization_time,
        observation_time, integrator, model, check_drift);
    // Note: scale_sig receives (evaluation + stabilization) time because both
    // segments are simulated before scoring starts.

    // Invalid scaling configuration (e.g. unsupported dt selection).
    if (!scaled_result_opt)
    {
        return std::nullopt;
    }

    const ScaledSigResult &scaled_result = *scaled_result_opt;

    // Post-synaptic neuron is advanced during candidate simulation.
    NeuronType model_neur = create_neur(false);
    using ChemicalSynapsisType = ChemicalSynapsis<NeuronType, NeuronType, Integrator, double>;

    typename ChemicalSynapsisType::ConstructorArgs syn_args{};
    if (use_i_slow && !use_i_fast)
    {
        // Disable fast branch physically in the synapse model.
        syn_args.params[ChemicalSynapsisType::gfast] = 0.0;
        syn_args.params[ChemicalSynapsisType::sfast] = 0.0;
        syn_args.params[ChemicalSynapsisType::Vfast] = 0.0;
    }
    else if (use_i_fast && !use_i_slow)
    {
        // Disable slow branch physically in the synapse model.
        syn_args.params[ChemicalSynapsisType::gslow] = 0.0;
        syn_args.params[ChemicalSynapsisType::Vslow] = 0.0;
        syn_args.params[ChemicalSynapsisType::k1] = 0.0;
        syn_args.params[ChemicalSynapsisType::k2] = 0.0;
        syn_args.params[ChemicalSynapsisType::sslow] = 0.0;
    }

    ChemicalSynapsisType synapse(create_neur(true), neur_v_var, model_neur, neur_v_var, syn_args, syn_model_step_factor);
    // Synapse takes (pre-neuron, pre-var, post-neuron, post-var, params, substeps).

    // First segment is stabilization-only; second segment is scored.
    const size_t stabilization_points = static_cast<size_t>(stabilization_time / csv_step);
    const size_t evaluation_points = scaled_result.sig.size() - stabilization_points;

    // Reused per-candidate output buffers.
    EvaluationISigBuffers buffers(evaluation_points, use_i_fast, use_i_slow);

    // Scored segment of the pre-synaptic voltage trace.
    const kfr::univector_ref<const double> evaluation_v_pre_sig = scaled_result.sig.slice(stabilization_points, evaluation_points);

    double v_post_min, v_post_max;
    if (model == NeuronModel::HINDMARSH_ROSE)
    {
        v_post_min = HindmarshRose::MIN;
        v_post_max = HindmarshRose::MAX;
    }
    else
    {
        throw std::invalid_argument("Unsupported neuron model");
    }

    const double evaluation_v_pre_min = kfr::minof(evaluation_v_pre_sig);
    const double evaluation_v_pre_max = kfr::maxof(evaluation_v_pre_sig);

    BOParamRanges ranges;
    // Decode bounds are built from observed signal + target current constraints.
    ranges.init(evaluation_v_pre_min, evaluation_v_pre_max, v_post_min, v_post_max,
                expected_i_min, expected_i_max,
                use_i_fast, use_i_slow,
                search_phase, fc);

    // Precompute reference scoring vectors/constants once.
    const ConstantEvaluationVals constant_evaluation_vals = calc_constant_evaluation_vals(
        evaluation_v_pre_sig,
        evaluation_v_pre_min,
        evaluation_v_pre_max,
        csv_step,
        fc,
        use_i_fast,
        use_i_slow,
        expected_i_min,
        expected_i_max,
        search_phase);

    // Dimensionality depends on enabled component subset.
    const int dim_in = static_cast<int>(count_params(use_i_fast, use_i_slow));
    const int iters = static_cast<int>(iterations);

    // Map user CLI configuration into Limbo dynamic parameters.
    Params::init_lhs::set_samples(static_cast<int>(initial_samples));
    Params::stop_maxiterations::set_iterations(iters);
    Params::bayes_opt_boptimizer::set_hp_period(std::max(BOPrivateConfig::HP_PERIOD_MIN, iters / BOPrivateConfig::HP_PERIOD_DIVISOR));
    Params::opt_rprop::set_iterations(BOPrivateConfig::RPROP_ITER_FACTOR * dim_in);
    Params::opt_nloptnograd::set_iterations(BOPrivateConfig::NLOPT_ITER_FACTOR * dim_in * dim_in);

    EvaluationFunctor<Integrator, NeuronType, ResetStateFuncType, GetVFuncType>::set_dim_in(dim_in);

    EvaluationFunctor<Integrator, NeuronType, ResetStateFuncType, GetVFuncType> functor(
        synapse,
        model_neur,
        scaled_result,
        buffers,
        constant_evaluation_vals,
        ranges,
        use_i_fast,
        use_i_slow,
        search_phase,
        stabilization_points,
        i_min,
        i_max,
        verbose,
        reset_state_neur,
        get_v_neur);

    BO_t opt;
    // Main optimization loop.
    opt.optimize(functor);

    if (verbose)
    {
        std::cout << "Best: " << opt.best_observation()(0) << std::endl;
    }
    // Decode best normalized sample to physical synapse parameters.
    const ChemicalSynapseParams best_candidate = decode_to_candidate(opt.best_sample(), ranges, use_i_fast, use_i_slow);

    return best_candidate;
}
