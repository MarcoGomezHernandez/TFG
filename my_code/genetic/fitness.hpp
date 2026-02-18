#ifndef FITNESS_H
#define FITNESS_H

#include <vector>
#include <array>
#include <ChemicalSynapsis.h>
#include "utils.hpp"
#include "scaling.hpp"

/*
 * Individual in the genetic algorithm population
 * Contains both parameters and fitness value
 */
struct Individual
{
    ChemicalSynapsisParams params;
    double fitness;
};

/*
 * Struct to hold precomputed signal statistics
 */
struct ConstantSignalFitnessVals
{
    std::vector<bool> model_up_states;
    std::vector<double> normalized_signal; // Normalized (and possibly inverted for antiphase) reference signal
    double bursts_seen;
    double norm_max_bursts_diff;
    double min_val;
    double max_val;
};

/*
 * Function to preprocess a signal and compute statistics
 * Returns ConstantSignalFitnessVals containing model_up_states, normalized_signal, bursts_seen, and norm_max_bursts_diff
 * search_phase: true for phase (normalized as-is), false for antiphase (normalized and y-inverted)
 */
ConstantSignalFitnessVals calc_const_signal_fitness_vals(const std::vector<double> &signal, double min_val, double max_val, bool search_phase);

/*
 * Template function to calculate fitnesses for multiple parameter sets
 * Simulates the neural model with given parameters and computes fitness against precomputed stats
 * Parameters: synapsis, neurons, individuals, scaled_result, stats1, search_phase, buffers, reset_state_neur, get_v_neur, start_index
 */
template <typename Integrator, typename NeuronType, size_t N, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    std::array<Individual, N> &individuals,
                    const ScaledSignalResult &scaled_result,
                    const ConstantSignalFitnessVals &stats1,
                    bool search_phase,
                    std::vector<double> &model_signal_buffer,
                    std::vector<double> &synapsis_signal_buffer,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur,
                    size_t start_index);

#include "fitness.tpp"

#endif
