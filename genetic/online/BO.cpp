#include <limbo/limbo.hpp>
#include "bidirectional_chemical_synapse_BO.h"

using namespace limbo;

namespace BOPrivateConfig
{
    // Padding length (in samples) is computed as: PAD_LEN_FACTOR * fs / fc.
    // Larger padding reduces filtfilt edge artifacts during scoring.
    static constexpr double PAD_LEN_FACTOR = 1.5;

    // Limbo hyper-parameter optimization (HP) period is clamped to at least this.
    static constexpr int HP_PERIOD_MIN = 1;
    // hp_period is computed as max(HP_PERIOD_MIN, iterations / HP_PERIOD_DIVISOR).
    static constexpr int HP_PERIOD_DIVISOR = 50;

    // Internal optimizer iteration multipliers (scaled by BO input dimension).
    static constexpr int RPROP_ITER_FACTOR = 30;
    static constexpr int NLOPT_ITER_FACTOR = 100;

    // Weights for the evaluation function of the BO (range vs shape).
    static constexpr double I_SHAPE_WEIGHT = 0.5;
    static constexpr double I_RANGE_WEIGH = 0.5;
}

namespace BOPrivateConstants
{
    static constexpr double TOTAL_WEIGHT = BOPrivateConfig::I_SHAPE_WEIGHT + BOPrivateConfig::I_RANGE_WEIGH;
}

struct Params
{
    struct bayes_opt_bobase : public defaults::bayes_opt_bobase
    {
        // Constrain the optimizer to [0, 1]^d (we decode to physical units manually).
        BO_PARAM(bool, bounded, true);
        // Enable Limbo statistics (samples, observations, acquisitions, ...).
        BO_PARAM(bool, stats_enabled, true);
    };

    struct bayes_opt_boptimizer : public defaults::bayes_opt_boptimizer
    {
        // Frequency (in iterations) for GP hyper-parameter updates.
        BO_DYN_PARAM(int, hp_period);
    };

    struct init_lhs : public defaults::init_lhs
    {
        // Number of initial Latin Hypercube samples.
        BO_DYN_PARAM(int, samples);
    };

    struct kernel : public defaults::kernel
    {
        // Observation noise assumed by the GP.
        BO_PARAM(double, noise, 0.05);
        // Allow Limbo to optimize noise as a GP hyper-parameter.
        BO_PARAM(bool, optimize_noise, true);
    };

    struct kernel_squared_exp_ard : public defaults::kernel_squared_exp_ard
    {
        // Initial dimension index used by ARD kernels (Limbo internal).
        BO_PARAM(int, k, 0);
        // Initial signal variance for the kernel.
        BO_PARAM(double, sigma_sq, 0.5);
    };

    struct acqui_gpucb : public defaults::acqui_gpucb
    {
        // Confidence parameter for GP-UCB acquisition.
        BO_PARAM(double, delta, 0.1);
    };

    struct opt_nloptgrad : public defaults::opt_nloptgrad
    {
        // Iterations for gradient-based acquisition optimization (scaled by dim).
        BO_DYN_PARAM(int, iterations);
        // Disable NLOpt stopping tolerances (let iteration count control runtime).
        BO_PARAM(double, fun_tolerance, -1);
        BO_PARAM(double, xrel_tolerance, -1);
    };

    struct opt_rprop : public defaults::opt_rprop
    {
        // Iterations for Rprop hyper-parameter optimization (scaled by dim).
        BO_DYN_PARAM(int, iterations);
        BO_PARAM(double, eps_stop, 1e-6);
    };

    struct stop_maxiterations : public defaults::stop_maxiterations
    {
        // Global BO iteration limit.
        BO_DYN_PARAM(int, iterations);
    };

    struct mean_constant : public defaults::mean_constant
    {
        // Constant mean function initial value.
        BO_PARAM(double, constant, 0.5);
    };
};

struct EvaluationFunctor
{
    // Bridges Limbo's evaluation API to the RTXI module.
    // Returns a 2D objective: [range_score, shape_score].
    EvaluationFunctor(BidirectionalChemicalSynapseBO &module_,
                      BOParamRanges &ranges_12_,
                      BOParamRanges &ranges_21_,
                      double fs_,
                      size_t effective_pad_12_,
                      size_t effective_pad_21_,
                      EvaluationPadBuffers &pad_buffers_,
                      double max_i_dist_12_,
                      double max_i_dist_21_,
                      size_t curr_synapse_idx_)
        : module(module_),
          ranges_12(ranges_12_),
          ranges_21(ranges_21_),
          fs(fs_),
          effective_pad_12(effective_pad_12_),
          effective_pad_21(effective_pad_21_),
          pad_buffers(pad_buffers_),
          max_i_dist_12(max_i_dist_12_),
          max_i_dist_21(max_i_dist_21_),
          curr_synapse_idx(curr_synapse_idx_)
    {
    }

    BidirectionalChemicalSynapseBO &module;
    BOParamRanges &ranges_12;
    BOParamRanges &ranges_21;
    double fs;
    size_t effective_pad_12, effective_pad_21;
    EvaluationPadBuffers &pad_buffers;
    double max_i_dist_12, max_i_dist_21;
    mutable size_t curr_synapse_idx;

    BO_DYN_PARAM(int, dim_in);
    BO_PARAM(int, dim_out, 1);

    Eigen::VectorXd operator()(const Eigen::VectorXd &x) const
    {
        // Decode normalized vector into physical synapse parameters.
        Candidate candidate = module.decode_to_candidate(x, ranges_12, ranges_21);

        // Evaluate by swapping params into RT, capturing signals, and scoring.
        ChemicalSynapseEvaluation evaluations = module.evaluate_candidate(candidate,
                                                                          fs,
                                                                          effective_pad_12,
                                                                          effective_pad_21,
                                                                          pad_buffers,
                                                                          max_i_dist_12,
                                                                          max_i_dist_21,
                                                                          curr_synapse_idx);

        const double y = ((evaluations.i_range_score * BOPrivateConfig::I_RANGE_WEIGH) + (evaluations.i_shape_score * BOPrivateConfig::I_SHAPE_WEIGHT)) / BOPrivateConstants::TOTAL_WEIGHT;
        return limbo::tools::make_vector(y);
    }
};

using Model_t = model::GP<
    Params,
    kernel::SquaredExpARD<Params>,
    mean::Constant<Params>,
    model::gp::KernelLFOpt<Params, opt::Rprop<Params>>>;

using Acqui_t = acqui::GP_UCB<Params, Model_t>;

using AcquiOpt_t = opt::NLOptGrad<Params, nlopt::LD_LBFGS>;

using Stat_t = boost::fusion::vector<
    stat::Samples<Params>,
    stat::Observations<Params>,
    stat::GPAcquisitions<Params>,
    stat::BestObservations<Params>,
    stat::AggregatedObservations<Params>>;

using Init_t = init::LHS<Params>;

using Stop_t = boost::fusion::vector<
    stop::MaxIterations<Params>,
    StopFunctor>;

using BO_t = bayes_opt::BOptimizer<
    Params,
    stopcrit<Stop_t>,
    modelfun<Model_t>,
    acquifun<Acqui_t>,
    statsfun<Stat_t>,
    initfun<Init_t>,
    acquiopt<AcquiOpt_t>>;

BO_DECLARE_DYN_PARAM(int, EvaluationFunctor, dim_in);
BO_DECLARE_DYN_PARAM(int, Params::bayes_opt_boptimizer, hp_period);
BO_DECLARE_DYN_PARAM(int, Params::init_lhs, samples);
BO_DECLARE_DYN_PARAM(int, Params::stop_maxiterations, iterations);
BO_DECLARE_DYN_PARAM(int, Params::opt_rprop, iterations);
BO_DECLARE_DYN_PARAM(int, Params::opt_nloptgrad, iterations);

static double decode_param(double x_val, const BOParamRanges::ParamRange &range)
{
    // Clamp to [0,1] and linearly map into [min, min+range].
    return range.min + (std::clamp(x_val, 0.0, 1.0) * range.range);
}

static void decode_to_params(const Eigen::VectorXd &x,
                             size_t &idx,
                             ChemicalSynapseParams &params,
                             unsigned int use_i_fast,
                             unsigned int use_i_slow,
                             const BOParamRanges &ranges)
{
    // Decode parameters in a fixed order; only enabled subsets consume dimensions.
    if (use_i_fast || use_i_slow)
    {
        params.e_syn = decode_param(x(idx++), ranges.e_syn);

        if (use_i_fast)
        {
            // Conductances are optimized in log-domain -> exp back to linear.
            params.g_fast = std::exp(decode_param(x(idx++), ranges.log_g_fast));
            params.s_fast = decode_param(x(idx++), ranges.s_fast);
            params.v_fast = decode_param(x(idx++), ranges.v_fast);
        }

        if (use_i_slow)
        {
            params.g_slow = std::exp(decode_param(x(idx++), ranges.log_g_slow));
            params.v_slow = decode_param(x(idx++), ranges.v_slow);
            params.k1 = std::exp(decode_param(x(idx++), ranges.log_k1));
            // R = k2/k1 is optimized in log-domain for stability.
            const double R = std::exp(decode_param(x(idx++), ranges.log_R));
            params.k2 = params.k1 * R;
            params.s_slow = decode_param(x(idx++), ranges.s_slow);
        }
    }
}

Candidate BidirectionalChemicalSynapseBO::decode_to_candidate(const Eigen::VectorXd &x,
                                                              const BOParamRanges &ranges_12,
                                                              const BOParamRanges &ranges_21)
{
    Candidate candidate;
    size_t idx = 0;
    decode_to_params(x, idx, candidate.params_12,
                     use_i_fast_12, use_i_slow_12, ranges_12);
    decode_to_params(x, idx, candidate.params_21,
                     use_i_fast_21, use_i_slow_21, ranges_21);
    return candidate;
}

static size_t count_params_one_direction(unsigned int use_i_fast,
                                         unsigned int use_i_slow)
{
    // Parameter count per direction:
    // - always: e_syn (1)
    // - fast: g_fast, s_fast, v_fast (3)
    // - slow: g_slow, v_slow, k1, R, s_slow (5)
    size_t num_params = 0;
    if (use_i_fast || use_i_slow)
    {
        num_params += 1;
        if (use_i_fast)
            num_params += 3;
        if (use_i_slow)
            num_params += 5;
    }
    return num_params;
}

void BidirectionalChemicalSynapseBO::NRT_BO(double period_t)
{
    // This function runs in the NRT BO thread. It must not call RT-only APIs.
    // It coordinates with the RT loop using atomics + a double-buffered params swap.
    const double fs = 1.0 / safe_divisor(period_t);           // period en ms, fs en KHz
    num_elements = static_cast<size_t>(evaluation_time * fs); // evaluation_time en ms
    if (num_elements < 1)
    {
        QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
        return;
    }
    const size_t pad_len_factor_fs = static_cast<size_t>(BOPrivateConfig::PAD_LEN_FACTOR * fs);

    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

    if (!use_syn_12 && !use_syn_21)
    {
        QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
        return;
    }

    size_t effective_pad_12 = 0;
    size_t padded_buff_size_12 = 0;
    if (use_syn_12)
    {
        // Effective padding scales inversely with fc (lower fc => longer impulse response).
        effective_pad_12 = std::min(num_elements - 1,
                                    static_cast<size_t>(pad_len_factor_fs / safe_divisor(fc_1)));
        padded_buff_size_12 = num_elements + (2 * effective_pad_12);
    }

    size_t effective_pad_21 = 0;
    size_t padded_buff_size_21 = 0;
    if (use_syn_21)
    {
        effective_pad_21 = std::min(num_elements - 1,
                                    static_cast<size_t>(pad_len_factor_fs / safe_divisor(fc_2)));
        padded_buff_size_21 = num_elements + (2 * effective_pad_21);
    }

    EvaluationPadBuffers pad_buffers(padded_buff_size_12, padded_buff_size_21);

    BOParamRanges ranges_12;
    if (use_syn_12)
    {
        // Build decoding ranges from observed voltages and target current bounds.
        ranges_12.init(v_min_1, v_max_1, v_min_2, v_max_2,
                       expected_i_min_12, expected_i_max_12,
                       use_i_fast_12, use_i_slow_12,
                       search_phase, fc_1);

        // Pre-allocate RT capture buffers (avoid allocation in execute()).
        v_sig_1.resize(num_elements);
        if (use_i_fast_12)
            i_fast_sig_12.resize(num_elements);
        if (use_i_slow_12)
            i_slow_sig_12.resize(num_elements);
    }

    BOParamRanges ranges_21;
    if (use_syn_21)
    {
        ranges_21.init(v_min_2, v_max_2, v_min_1, v_max_1,
                       expected_i_min_21, expected_i_max_21,
                       use_i_fast_21, use_i_slow_21,
                       search_phase, fc_2);
        v_sig_2.resize(num_elements);
        if (use_i_fast_21)
            i_fast_sig_21.resize(num_elements);
        if (use_i_slow_21)
            i_slow_sig_21.resize(num_elements);
    }

    const int dim_in = static_cast<int>(count_params_one_direction(use_i_fast_12, use_i_slow_12) + count_params_one_direction(use_i_fast_21, use_i_slow_21));
    const int iters = static_cast<int>(iterations);

    // Map module GUI config into Limbo dynamic parameters.
    Params::init_lhs::set_samples(static_cast<int>(initial_samples));
    Params::stop_maxiterations::set_iterations(iters);
    Params::bayes_opt_boptimizer::set_hp_period(std::max(BOPrivateConfig::HP_PERIOD_MIN, iters / BOPrivateConfig::HP_PERIOD_DIVISOR));
    Params::opt_rprop::set_iterations(BOPrivateConfig::RPROP_ITER_FACTOR * dim_in);
    Params::opt_nloptgrad::set_iterations(BOPrivateConfig::NLOPT_ITER_FACTOR * dim_in);

    EvaluationFunctor::set_dim_in(dim_in);

    double max_i_dist_12 = 0.0;
    if (use_syn_12)
    {
        // If only one component is enabled, use the user-provided expected range directly.
        // If both are enabled, evaluate_sigs_one_direction computes per-component ranges.
        if (use_i_fast_12 != use_i_slow_12)
        {
            max_i_dist_12 = calculate_expected_i_max_dist(expected_i_min_12, expected_i_max_12);
        }
    }
    double max_i_dist_21 = 0.0;
    if (use_syn_21)
    {
        if (use_i_fast_21 != use_i_slow_21)
        {
            max_i_dist_21 = calculate_expected_i_max_dist(expected_i_min_21, expected_i_max_21);
        }
    }

    EvaluationFunctor functor(*this,
                              ranges_12,
                              ranges_21,
                              fs,
                              effective_pad_12,
                              effective_pad_21,
                              pad_buffers,
                              max_i_dist_12,
                              max_i_dist_21,
                              // Start from the current RT-visible params slot.
                              synapse_idx.load(std::memory_order_relaxed));

    QMetaObject::invokeMethod(this, "set_evaluations_completed", Qt::QueuedConnection,
                              Q_ARG(double, 0));

    BO_t opt;
    try
    {
        // Main optimization loop (stops when iterations are reached or stop flag triggers).
        opt.optimize(functor);
    }
    catch (const StopEvaluation &)
    {
        QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
        return;
    }

    if (stop_BO.load(std::memory_order_relaxed))
    {
        QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
        return;
    }

    Candidate best_candidate = decode_to_candidate(opt.best_sample(),
                                                   ranges_12,
                                                   ranges_21);

    const size_t curr_synapse_idx = functor.curr_synapse_idx;
    // Ensure RT has observed the current slot before swapping to the other one.
    if (!wait_until_RT_read_idx_or_stop(curr_synapse_idx))
    {
        QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
        return;
    }

    const size_t new_synapse_idx = 1 - curr_synapse_idx;
    // Write the best candidate into the inactive slot, then publish via release store.
    copy_selected_synapse_params(params_12[new_synapse_idx],
                                 best_candidate.params_12,
                                 use_i_fast_12,
                                 use_i_slow_12);
    copy_selected_synapse_params(params_21[new_synapse_idx],
                                 best_candidate.params_21,
                                 use_i_fast_21,
                                 use_i_slow_21);
    synapse_idx.store(new_synapse_idx, std::memory_order_release);

    QMetaObject::invokeMethod(this, "stop_BO_event_async", Qt::QueuedConnection);
}
