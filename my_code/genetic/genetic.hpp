#ifndef GENETIC_H
#define GENETIC_H

/*
 * Public genetic algorithm configuration constants
 */
namespace GeneticPublicConfig
{
    // Fixed values for gfast and gslow
    inline constexpr double GFAST_FIXED = 0.046;
    inline constexpr double GSLOW_FIXED = 0.208;
}

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
                   CreateFuncType create_neuron,
                   ResetStateFuncType reset_state_neur,
                   GetVFuncType get_v_neur);

#include "genetic.tpp"

#endif // GENETIC_H
