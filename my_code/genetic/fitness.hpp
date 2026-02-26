#ifndef FITNESS_H
#define FITNESS_H

#include <array>
#include <kfr/all.hpp>
#include <ChemicalSynapsis.h>
#include "utils.hpp"
#include "scaling.hpp"

struct Individual
{
    ChemicalSynapsisParams params;
    double fitness;
};

struct ConstantSignalFitnessVals
{
    kfr::univector<double> normalized_signal_to_fit;
    kfr::univector<double> smoothed_signal_to_fit;
    kfr::univector<double> norm_signal_to_fit_ifast;
    kfr::univector<double> norm_signal_to_fit_islow;
    double max_v_comp_distance;
};

struct SignalBuffers
{
    kfr::univector<double> model_signal;
    kfr::univector<double> synapsis_signal;
    kfr::univector<double> ifast_signal;
    kfr::univector<double> islow_signal;

    kfr::univector<double> kfr_padded;
};

ConstantSignalFitnessVals calc_const_signal_fitness_vals(const kfr::univector<double> &signal, double min_val, double max_val, bool search_phase, size_t avg_smooth_points, size_t start_index, double fs, double fc, SignalBuffers &buffers, bool use_ifast, bool use_islow);

template <typename Integrator, typename NeuronType, size_t N, ResetStateFunc<NeuronType> ResetStateFuncType, GetVFunc<NeuronType> GetVFuncType>
void calc_fitnesses(ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> &synapsis,
                    NeuronType &model_neur,
                    std::array<Individual, N> &individuals,
                    const ScaledSignalResult &scaled_result,
                    const ConstantSignalFitnessVals &stats1,
                    bool search_phase,
                    SignalBuffers &buffers,
                    ResetStateFuncType reset_state_neur,
                    GetVFuncType get_v_neur,
                    size_t ind_start_index,
                    size_t signal_start_index,
                    size_t avg_smooth_points_model,
                    bool use_ifast,
                    bool use_islow);

#include "fitness.tpp"

#endif
