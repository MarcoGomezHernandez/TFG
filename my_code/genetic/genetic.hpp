#ifndef GENETIC_H
#define GENETIC_H

#include <string>
#include "utils.hpp"
#include "fitness.hpp"

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
Individual genetic(const std::string &csv_path,
                   size_t column_index,
                   double csv_step,
                   double start_time,
                   double use_time,
                   NumericIntegrator integrator,
                   NeuronModel model,
                   bool search_phase,
                   bool check_drift,
                   CreateFuncType create_neur,
                   ResetStateFuncType reset_state_neur,
                   GetVFuncType get_v_neur,
                   typename NeuronType::VarType neur_state_var,
                   int syn_model_step_factor);

#include "genetic.tpp"

#endif // GENETIC_H
