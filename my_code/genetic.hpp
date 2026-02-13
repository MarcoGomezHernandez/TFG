#ifndef GENETIC_H
#define GENETIC_H

#include <string>
#include <array>
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <utility>
#include <stdexcept>

#include "utils.hpp"
#include "scaling.hpp"
#include "fitness.hpp"
#include <ChemicalSynapsis.h>

struct BestResult
{
    double fitness;
    ChemicalSynapsisParams params;
};

/*
 * Main genetic algorithm function template.
 * NeuronType: the neuron wrapper type
 * Integrator: the numeric integrator type
 * CreateFuncType: callable returning NeuronType
 * ResetStateFuncType: callable taking NeuronType& and resetting state
 * GetVFuncType: callable taking const NeuronType& and returning double (voltage)
 */
template <typename Integrator, typename NeuronType,
          CreateFunc<NeuronType> CreateFuncType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
BestResult genetic(const std::string &csv_path,
                   size_t column_index,
                   double csv_step,
                   double start_time,
                   double use_time,
                   NumericIntegrator integrator_enum,
                   NeuronModel model_enum,
                   bool search_phase,
                   bool check_drift,
                   CreateFuncType create_neuron,
                   ResetStateFuncType reset_state_neur,
                   GetVFuncType get_v_neur);

#endif // GENETIC_H
