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
    if (argc < 9)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <csv_path> <column_index> <csv_step> <start_time> <use_time> <search_phase> <check_drift> <steps>"
                  << std::endl;
        return 1;
    }

    const std::string csv_path = argv[1];
    const size_t column_index = static_cast<size_t>(std::atoi(argv[2]));
    const double csv_step = std::atof(argv[3]);
    const double start_time = std::atof(argv[4]);
    const double use_time = std::atof(argv[5]);
    const bool search_phase = (std::atoi(argv[6]) != 0);
    const bool check_drift = (std::atoi(argv[7]) != 0);
    const int steps = std::atoi(argv[8]);

    try
    {
        const Individual &best = genetic<Integrator, NeuronType>(
            csv_path, column_index, csv_step, start_time, use_time,
            NumericIntegrator::RK4,
            NeuronModel::HINDMARSH_ROSE,
            search_phase,
            check_drift,
            create_hindmarsh_rose<Integrator>,
            reset_state_hindmarsh_rose<Integrator>,
            get_v_hindmarsh_rose<Integrator>,
            NeuronType::x,
            steps);

        const ChemicalSynapsisParams &params = best.params;

        std::cout << "=== Best Individual ===" << std::endl;
        std::cout << "Fitness: " << best.fitness << std::endl;
        std::cout << "gfast:  " << GeneticPublicConfig::GFAST_FIXED << " (fixed)" << std::endl;
        std::cout << "Esyn:   " << params.Esyn << std::endl;
        std::cout << "sfast:  " << params.sfast << std::endl;
        std::cout << "Vfast:  " << params.Vfast << std::endl;
        std::cout << "Vslow:  " << params.Vslow << std::endl;
        std::cout << "gslow:  " << GeneticPublicConfig::GSLOW_FIXED << " (fixed)" << std::endl;
        std::cout << "k1:     " << params.k1 << std::endl;
        std::cout << "k2:     " << params.k2 << std::endl;
        std::cout << "sslow:  " << params.sslow << std::endl;
    }
    catch (const std::runtime_error &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
