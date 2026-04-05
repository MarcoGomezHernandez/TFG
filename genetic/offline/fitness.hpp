#ifndef FITNESS_H
#define FITNESS_H

#include <array>
#include <span>
#include <kfr/all.hpp>
using namespace kfr;
#include <ChemicalSynapsis.h>
#include "utils.hpp"
#include "scaling.hpp"

namespace FitnessPublicConfig
{
  inline constexpr double EXPECTED_I_MIN_ANTIPHASE = 0.2;
  inline constexpr double EXPECTED_I_MAX_ANTIPHASE = 0.8;

  inline constexpr double EXPECTED_I_MIN_PHASE = -1.4;
  inline constexpr double EXPECTED_I_MAX_PHASE = -0.4;
}

struct Individual
{
  ChemicalSynapsisVariationParams params;
  double fitness;
};

struct ConstantSigFitnessVals
{
  univector<double> ifast_sig_centered_to_fit;
  double ifast_sig_factor_to_fit;
  univector<double> islow_sig_centered_to_fit;
  double islow_sig_factor_to_fit;

  double ifast_sig_min_to_fit;
  double ifast_sig_max_to_fit;
  double ifast_sig_range_to_fit;

  double islow_sig_min_to_fit;
  double islow_sig_max_to_fit;
  double islow_sig_range_to_fit;
};

struct SigBuffers
{
  SigBuffers(size_t size_to_reserve,
             bool use_ifast,
             bool use_islow)
  {
    if (use_ifast)
      ifast_sig.resize(size_to_reserve);
    if (use_islow)
      islow_sig.resize(size_to_reserve);
  }

  univector<double> ifast_sig;
  univector<double> islow_sig;
};

ConstantSigFitnessVals calc_const_sig_fitness_vals(const univector_ref<const double> &vpre_sig,
                                                   double csv_step,
                                                   bool use_ifast,
                                                   bool use_islow,
                                                   bool search_phase);

template <typename Integrator, typename NeuronType, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    std::span<Individual> individuals,
                    const ScaledSigResult &scaled_result,
                    const ConstantSigFitnessVals &const_vpre_sig_fitness_vals,
                    bool search_phase,
                    SigBuffers &buffers,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur,
                    size_t vpre_sig_start_i,
                    bool use_ifast,
                    bool use_islow);

#include "fitness.tpp"

#endif
