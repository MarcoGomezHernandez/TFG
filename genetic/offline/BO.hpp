#ifndef BO_H
#define BO_H

#include <limbo/limbo.hpp>
#include <Eigen/Core>
#include <cstddef>
#include <optional>
#include <string>
#include "utils.hpp"
#include "evaluation.hpp"

// Entry point for offline Bayesian Optimization of one chemical synapse direction.
//
// Workflow summary:
// 1) Read and scale a reference pre-synaptic voltage trace from CSV.
// 2) Build BO parameter ranges from signal/model constraints.
// 3) Run Limbo BO to maximize a weighted score (range + shape).
// 4) Decode and return the best physical synapse parameters.
template <typename Integrator, typename NeuronType,
          CreateFunc<NeuronType> CreateFuncType,
          ResetStateFunc<NeuronType> ResetStateFuncType,
          GetVFunc<NeuronType> GetVFuncType>
std::optional<ChemicalSynapseParams> BO(const std::string &csv_path,
                                        size_t column_idx,
                                        double csv_step,
                                        double start_time,
                                        double stabilization_time,
                                        double evaluation_time,
                                        double observation_time,
                                        size_t initial_samples,
                                        size_t iterations,
                                        NumericIntegrator integrator,
                                        NeuronModel model,
                                        bool search_phase,
                                        bool check_drift,
                                        SynComponent syn_component,
                                        CreateFuncType create_neur,
                                        ResetStateFuncType reset_state_neur,
                                        GetVFuncType get_v_neur,
                                        typename NeuronType::VarType neur_v_var,
                                        int syn_model_step_factor,
                                        double fc,
                                        double expected_i_min,
                                        double expected_i_max,
                                        double i_min,
                                        double i_max,
                                        bool verbose);

#include "BO.tpp"

#endif
