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
    std::vector<double> normalized_signal_to_fit; // Normalized (and possibly flipped for antiphase) reference signal
    std::vector<double> smoothed_signal_to_fit;   // Smoothed signal (possibly flipped for antiphase)
    double max_phase_distance;              // Precomputed: signal_size * (max_smoothed_val - min_smoothed_val)
};

/*
 * Compute fitness score for a single individual's synapsis signal against precomputed reference signal
 * Combines phase similarity and biological similarity into a weighted score
 * Parameters: signal, min_val, max_val, search_phase, avg_smooth_points
 */
ConstantSignalFitnessVals calc_const_signal_fitness_vals(const std::vector<double> &signal, double min_val, double max_val, bool search_phase, size_t avg_smooth_points);

/*
 * Template function to calculate fitnesses for multiple parameter sets
 * Simulates the neural model with given parameters and computes fitness against precomputed stats
 * Parameters: synapsis, neurons, individuals, scaled_result, stats1, search_phase, buffers, reset_state_neur, get_v_neur, start_index, avg_smooth_points_model
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
                    size_t start_index,
                    size_t avg_smooth_points_model);

#include "fitness.tpp"

#endif
