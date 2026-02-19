#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>

#include "genetic.hpp"
#include "utils.hpp"
#include "scaling.hpp"

#include <RungeKutta4.h>
#include <DifferentialNeuronWrapper.h>
#include <SystemWrapper.h>
#include <HindmarshRoseModel.h>
#include <ChemicalSynapsis.h>

typedef RungeKutta4 Integrator;
typedef HindmarshRoseNeuron<Integrator> NeuronType;

int main(int argc, char *argv[])
{
    if (argc < 10)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <csv_path> <column_index> <csv_step> <start_time> <use_time> <search_phase> <check_drift> <syn_model_step_factor> <syn_component>"
                  << std::endl;
        std::cerr << "  syn_component: 0=ifast  1=islow  2=both" << std::endl;
        return 1;
    }

    const std::string csv_path = argv[1];
    const size_t column_index = static_cast<size_t>(std::atoi(argv[2]));
    const double csv_step = std::atof(argv[3]);
    const double start_time = std::atof(argv[4]);
    const double use_time = std::atof(argv[5]);
    const bool search_phase = (std::atoi(argv[6]) == 1);
    const bool check_drift = (std::atoi(argv[7]) == 1);
    const int syn_model_step_factor = std::atoi(argv[8]);
    const SynComponent syn_component = static_cast<SynComponent>(std::atoi(argv[9]));

    try
    {
        const Individual &best = genetic<Integrator, NeuronType>(
            csv_path, column_index, csv_step, start_time, use_time,
            NumericIntegrator::RK4,
            NeuronModel::HINDMARSH_ROSE,
            search_phase,
            check_drift,
            syn_component,
            create_hindmarsh_rose<Integrator>,
            reset_state_hindmarsh_rose<Integrator>,
            get_v_hindmarsh_rose<Integrator>,
            NeuronType::x,
            syn_model_step_factor,
            true);

        // Output moved to `genetic()` when verbose is enabled.
    }
    catch (const std::runtime_error &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
