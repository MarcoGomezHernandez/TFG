#ifndef FITNESS_H
#define FITNESS_H

#include <vector>
#include <array>
#include "utils.h"
#include "scaling.h"
#include <ChemicalSynapsis.h>

/*
 * Struct to hold precomputed signal statistics
 */
struct ConstantSignalFitnessVals
{
    std::vector<bool> up_states;
    double bursts_seen;
    double norm_max_bursts_diff;
    double min_val;
    double max_val;
};

/*
 * Function to preprocess a signal and compute statistics
 * Returns ConstantSignalFitnessVals containing up_states, bursts_seen, and norm_max_bursts_diff
 */
ConstantSignalFitnessVals calc_const_signal_vals(const std::vector<double> &signal, double min_val, double max_val);

/*
 * Template function to calculate fitnesses for multiple parameter sets
 * Simulates the neural model with given parameters and computes fitness against precomputed stats
 * Parameters: synapsis, neurons, params_individuals, scaled_result, initial_state, stats1, search_phase, buffers, reset_state_neur, get_v_neur, set_v_neur
 */
template <typename Integrator, typename NeuronType, size_t N, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    const std::array<ChemicalSynapsisParams, N> &params_individuals,
                    const ScaledSignalResult &scaled_result,
                    const ConstantSignalFitnessVals &stats1,
                    bool search_phase,
                    std::vector<double> &model_signal_buffer,
                    std::array<double, N> &fitnesses_buffer,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur);

#endif
