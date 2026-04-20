#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <DifferentialNeuronWrapper.h>
#include <HindmarshRoseModel.h>
#include <RungeKutta4.h>
#include <SystemWrapper.h>
#include "aux/utils.hpp"

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<HindmarshRoseModel<double>>, Integrator> HRNeuron;

int main()
{
    HRNeuron n = create_hindmarsh_rose<Integrator>(false);
    reset_state_hindmarsh_rose(n);

    const double step = 0.001; // Paso de integración (estabilidad numérica)
    const int subsample = 100; // Guarda 1 de cada 100 iteraciones

    const double stabilization_time = 500.0;
    const double simulation_time = 7000.0;

    // 1. Estabilización (Burn-in transitorio)
    for (double t = 0; t < stabilization_time; t += step)
    {
        n.step(step);
    }

    std::ofstream out("data/hindmarsh-rose_modified.csv");
    out << "t,x,y,z\n";

    // 2. Simulación con subsampling
    int step_count = 0;
    for (double t = 0; t < simulation_time; t += step)
    {
        if (step_count % subsample == 0)
        {
            // El tiempo registrado cuenta desde que termina la estabilización
            out << t << "," << n.get(HRNeuron::x) << ","
                << n.get(HRNeuron::y) << "," << n.get(HRNeuron::z) << "\n";
        }
        n.step(step);
        step_count++;
    }

    out.close();
    return 0;
}
