#include <iostream>
#include <string>
#include <cstdlib>
#include <exception>
#include <optional>
#include <chrono>

#include "BO.hpp"
#include "utils.hpp"
#include "scaling.hpp"

typedef RungeKutta4 Integrator;
typedef HindmarshRoseNeuron<Integrator> NeuronType;

int main(int argc, char *argv[])
{

    if (argc < 20)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <csv_path> <column_idx> <csv_step (ms)> <start_time (ms)> <stabilization_time (ms)> <evaluation_time (ms)> <observation_time (ms)> <initial_samples> <iterations> <search_phase> <check_drift> <syn_model_step_factor> <syn_component> <cutoff_frequency (kHz)> <expected_i_min> <expected_i_max> <i_min> <i_max> <verbose>"
                  << std::endl;
        std::cerr << "  If we consider the post neuron in V, the currents are in nA" << std::endl;
        std::cerr << "  syn_component: 0=ifast  1=islow  2=both" << std::endl;
        return 1;
    }

    // Parsea los argumentos de línea de comandos
    const std::string csv_path = argv[1];
    const size_t column_idx = static_cast<size_t>(std::atoi(argv[2]));
    const double csv_step = std::atof(argv[3]);                                        // Paso temporal del CSV (ms)
    const double start_time = std::atof(argv[4]);                                      // Tiempo de inicio de lectura del CSV (ms)
    const double stabilization_time = std::atof(argv[5]);                              // Tiempo de estabilización antes de evaluar (ms)
    const double evaluation_time = std::atof(argv[6]);                                 // Tiempo de evaluación de la sinapsis (ms)
    const double observation_time = std::atof(argv[7]);                                // Tiempo para calcular estadísticas de la señal (ms)
    const size_t initial_samples = static_cast<size_t>(std::atoi(argv[8]));            // Muestras iniciales LHS
    const size_t iterations = static_cast<size_t>(std::atoi(argv[9]));                 // Iteraciones de BO
    const bool search_phase = (std::atoi(argv[10]) == 1);                              // 1 = fase (excitadora), 0 = antifase (inhibidora)
    const bool check_drift = (std::atoi(argv[11]) == 1);                               // Corrección de drift en la señal
    const int syn_model_step_factor = std::atoi(argv[12]);                             // Sub-pasos del modelo sináptico por paso
    const SynComponent syn_component = static_cast<SynComponent>(std::atoi(argv[13])); // 0=ifast, 1=islow, 2=both
    const double fc = std::atof(argv[14]);                                             // Frecuencia de corte Butterworth (kHz)
    const double expected_i_min = std::atof(argv[15]);                                 // Corriente mínima esperada (nA)
    const double expected_i_max = std::atof(argv[16]);                                 // Corriente máxima esperada (nA)
    const double i_min = std::atof(argv[17]);                                          // Clamp mínimo de salida (nA)
    const double i_max = std::atof(argv[18]);                                          // Clamp máximo de salida (nA)
    const bool verbose = (std::atoi(argv[19]) == 1);

    const std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();

    try
    {
        // Lanza la optimización bayesiana offline
        const std::optional<ChemicalSynapseParams> best_params_opt = BO<Integrator, NeuronType>(
            csv_path, column_idx, csv_step, start_time, stabilization_time, evaluation_time,
            observation_time,
            initial_samples,
            iterations,
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
            fc,
            expected_i_min,
            expected_i_max,
            i_min,
            i_max,
            verbose);

        const std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = t_end - t_start;

        if (!best_params_opt)
        {

            std::cerr << "Error: BO failed to produce parameters" << std::endl;
            std::cerr << "BO execution time (until failure): " << elapsed.count() << " s" << std::endl;
            return 1;
        }

        // Imprime los mejores parámetros encontrados
        const ChemicalSynapseParams &best_params = *best_params_opt;

        std::cout << "    Esyn_values = [" << best_params.e_syn << "]" << std::endl;
        std::cout << "    gfast_values = [" << best_params.g_fast << "]" << std::endl;
        std::cout << "    sfast_values = [" << best_params.s_fast << "]" << std::endl;
        std::cout << "    Vfast_values = [" << best_params.v_fast << "]" << std::endl;
        std::cout << "    gslow_values = [" << best_params.g_slow << "]" << std::endl;
        std::cout << "    k1_values = [" << best_params.k1 << "]" << std::endl;
        std::cout << "    k2_values = [" << best_params.k2 << "]" << std::endl;
        std::cout << "    sslow_values = [" << best_params.s_slow << "]" << std::endl;
        std::cout << "    Vslow_values = [" << best_params.v_slow << "]" << std::endl;
        std::cout << "BO execution time: " << elapsed.count() << " s" << std::endl;
    }
    catch (const std::exception &e)
    {

        const std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = t_end - t_start;
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "BO execution time (until error): " << elapsed.count() << " s" << std::endl;
        return 1;
    }

    return 0;
}
