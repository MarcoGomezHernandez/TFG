#include <iostream>
#include <string>
#include <cstdlib>

#include "genetic.hpp"
#include "utils.hpp"
#include "scaling.hpp"

#include <RungeKutta4.h>
#include <DifferentialNeuronWrapper.h>
#include <SystemWrapper.h>
#include <HindmarshRoseModel.h>
#include <ChemicalSynapsis.h>

int main(int argc, char *argv[])
{
    if (argc < 6)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <csv_path> <column_index> <csv_step> <start_time> <use_time> [search_phase=1] [check_drift=0]"
                  << std::endl;
        return 1;
    }

    const std::string csv_path = argv[1];
    const size_t column_index = static_cast<size_t>(std::atoi(argv[2]));
    const double csv_step = std::atof(argv[3]);
    const double start_time = std::atof(argv[4]);
    const double use_time = std::atof(argv[5]);
    const bool search_phase = (argc > 6) ? (std::atoi(argv[6]) != 0) : true;
    const bool check_drift = (argc > 7) ? (std::atoi(argv[7]) != 0) : false;

    typedef RungeKutta4 Integrator;
    typedef HindmarshRoseNeuron<Integrator> NeuronType;

    genetic<Integrator, NeuronType>(
        csv_path, column_index, csv_step, start_time, use_time,
        NumericIntegrator::RK4,
        NeuronModel::HINDMARSH_ROSE,
        search_phase,
        check_drift,
        create_hindmarsh_rose<Integrator>,
        reset_state_hindmarsh_rose<Integrator>,
        get_v_hindmarsh_rose<Integrator>);

    return 0;
}
