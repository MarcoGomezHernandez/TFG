#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <DifferentialNeuronWrapper.h>
#include <HindmarshRoseModel.h>
#include <RungeKutta4.h>
#include <SystemWrapper.h>

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<HindmarshRoseModel<double>>, Integrator> HRNeuron;

int main()
{
    HRNeuron::ConstructorArgs args;

    args.params[HRNeuron::e] = 3.281;
    args.params[HRNeuron::mu] = 0.0021;
    args.params[HRNeuron::S] = 1.0;
    args.params[HRNeuron::a] = 1.0;
    args.params[HRNeuron::b] = 3.0;
    args.params[HRNeuron::c] = 1.0;
    args.params[HRNeuron::d] = 5.0;
    args.params[HRNeuron::xr] = -1.6;
    args.params[HRNeuron::vh] = 0.1;

    HRNeuron n(args);

    n.set(HRNeuron::x, -0.712841);
    n.set(HRNeuron::y, -1.93688);
    n.set(HRNeuron::z, 3.16568);

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
