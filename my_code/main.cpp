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

    try
    {
        Individual best = genetic<Integrator, NeuronType>(
            csv_path, column_index, csv_step, start_time, use_time,
            NumericIntegrator::RK4,
            NeuronModel::HINDMARSH_ROSE,
            search_phase,
            check_drift,
            create_hindmarsh_rose<Integrator>,
            reset_state_hindmarsh_rose<Integrator>,
            get_v_hindmarsh_rose<Integrator>);

        std::cout << "=== Best Individual ===" << std::endl;
        std::cout << "Fitness: " << best.fitness << std::endl;
        std::cout << "gfast:  " << best.params.gfast << std::endl;
        std::cout << "Esyn:   " << best.params.Esyn << std::endl;
        std::cout << "sfast:  " << best.params.sfast << std::endl;
        std::cout << "Vfast:  " << best.params.Vfast << std::endl;
        std::cout << "Vslow:  " << best.params.Vslow << std::endl;
        std::cout << "gslow:  " << best.params.gslow << std::endl;
        std::cout << "k1:     " << best.params.k1 << std::endl;
        std::cout << "k2:     " << best.params.k2 << std::endl;
        std::cout << "sslow:  " << best.params.sslow << std::endl;
    }
    catch (const std::runtime_error &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
