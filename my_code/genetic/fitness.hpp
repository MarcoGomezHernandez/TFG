#ifndef FITNESS_H
#define FITNESS_H

#include <array>
#include <kfr/all.hpp>
using namespace kfr;
#include <ChemicalSynapsis.h>
#include "utils.hpp"
#include "scaling.hpp"

struct Individual
{
  ChemicalSynapsisParams params;
  double fitness;
};

struct ConstantSigFitnessVals
{
  univector<double> i_sig_centered_to_fit;
  double i_sig_stddev_to_fit;
  univector<double> ifast_sig_centered_to_fit;
  double ifast_sig_stddev_to_fit;
  univector<double> islow_sig_centered_to_fit;
  double islow_sig_stddev_to_fit;
  univector<double> smoothed_vpre_sig_to_fit;
  double max_v_comp_distance;
};

struct SigBuffers
{
  univector<double> vpost_sig;
  univector<double> smoothed_vpost_sig;
  univector<double> i_sig;
  univector<double> ifast_sig;
  univector<double> islow_sig;
};

ConstantSigFitnessVals calc_const_sig_fitness_vals(const univector<double> &vpre_sig,
                                                   double min,
                                                   double max,
                                                   bool search_phase,
                                                   size_t avg_smooth_points,
                                                   double pts_burst_real,
                                                   bool use_ifast,
                                                   bool use_islow);

template <typename Integrator, typename NeuronType, size_t N, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    std::array<Individual, N> &individuals,
                    const ScaledSigResult &scaled_result,
                    const ConstantSigFitnessVals &const_vpre_sig_fitness_vals,
                    bool search_phase,
                    SigBuffers &buffers,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur,
                    size_t ind_start_i,
                    size_t vpre_sig_start_i,
                    size_t avg_smooth_points,
                    bool use_ifast,
                    bool use_islow);

#include "fitness.tpp"

#endif
