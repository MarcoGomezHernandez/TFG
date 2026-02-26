#ifndef GENETIC_H
#define GENETIC_H

#include <string>
#include "utils.hpp"
#include "fitness.hpp"

template <typename Integrator, typename NeuronType,
          CreateFunc<NeuronType> CreateFuncType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
Individual genetic(const std::string &csv_path,
                   size_t column_i,
                   double csv_step,
                   double start_time,
                   double use_time,
                   double stabilization_time,
                   NumericIntegrator integrator,
                   NeuronModel model,
                   bool search_phase,
                   bool check_drift,
                   SynComponent syn_component,
                   CreateFuncType create_neur,
                   ResetStateFuncType reset_state_neur,
                   GetVFuncType get_v_neur,
                   typename NeuronType::VarType neur_state_var,
                   int syn_model_step_factor,
                   bool verbose);

#include "genetic.tpp"

#endif
