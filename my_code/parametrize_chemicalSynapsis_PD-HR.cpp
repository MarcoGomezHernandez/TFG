/*************************************************************

Copyright (c) 2025,
Alicia Garrido Peña <alicia.garrido@uam.es>,
Pablo Sanchez Martin <pablo.sanchezm@uam.es>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
    * Neither the name of the author nor the names of his contributors
      may be used to endorse or promote products derived from this
      software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*************************************************************/

#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <filesystem>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include <DifferentialNeuronWrapper.h>
#include <ChemicalSynapsis.h>
#include <HindmarshRoseModel.h>
#include <SystemWrapper.h>
#include <RungeKutta4.h>
#include "scaling.h"

using namespace std;

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<HindmarshRoseModel<double>>, Integrator> HR;
typedef ChemicalSynapsis<HR, HR, Integrator, double> Synapsis;

// Parse YAML section into a map
std::map<std::string,double> parse_file(const std::string& file, const std::string& section)
{
    std::map<std::string,double> config_values;
    YAML::Node config = YAML::LoadFile(file);
    if (!config[section]) {
        std::cerr << "Section '" << section << "' not found in " << file << std::endl;
        return config_values;
    }
    for (auto it = config[section].begin(); it != config[section].end(); ++it) {
        std::string key = it->first.as<std::string>();
        double value = it->second.as<double>();
        config_values[key] = value;
    }
    return config_values;
}

int main(int argc, char **argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <yaml_file> <output_file> <simulation_time> <csv_file> [step]\n";
        return 1;
    }

    // --- Command line arguments ---
    std::string yaml_file   = argv[1];
    std::string output_file = argv[2];
    double simulation_time  = std::atof(argv[3]);
    double real_step = 0.0;
    std::string csv_file;
    if (argc == 6) {
        real_step = std::atof(argv[4]);
        csv_file = argv[5];
    } else {
        csv_file = argv[4];
    }

    cout << "YAML file: " << yaml_file << "\n";
    cout << "Output file: " << output_file << "\n";
    cout << "Simulation time: " << simulation_time << "\n";
    cout << "CSV file: " << csv_file << "\n";
    cout << "Real step: " << real_step << "\n";

    // --- Initialize postsynaptic neuron HR1 ---
    HR::ConstructorArgs args;
    auto param_names = HR::ParamNames();
    auto v_args = parse_file(yaml_file, "HR1");
    for(int i=0; i < HR::n_parameters; i++)
        args.params[i] = v_args[param_names[i]];

    HR h1(args);
    h1.set(HR::x, -0.712841);
    h1.set(HR::y, -1.93688);
    h1.set(HR::z, 3.16568);

    // --- Initialize dummy presynaptic neuron for CSV ---
    HR::ConstructorArgs args_csv;
    HR HR_csv(args_csv);
    HR_csv.set(HR::x, 0.0); // initial value

    // --- Initialize synapse (CSV → HR1) ---
    Synapsis::ConstructorArgs syn_args;
    v_args = parse_file(yaml_file, "Chemical-HR1-HR2");
    auto syn_param_names = Synapsis::ParamNames();
    for(int i=0; i < Synapsis::n_parameters; i++)
        syn_args.params[i] = v_args[syn_param_names[i]];

    Synapsis s(HR_csv, HR::x, h1, HR::x, syn_args, 1);

    // --- Open output file ---
    ofstream out(output_file);
    if (!out.is_open()) {
        cerr << "Error: could not open " << output_file << " for writing.\n";
        return 1;
    }
    out << "Time Vpre Vpost i ifast islow\n";

    // --- Distance file ---
    string dist_file = "distances.csv";
    bool exists = filesystem::exists(dist_file);
    ofstream dist_out(dist_file, ios::app);
    if (!dist_out.is_open()) {
        cerr << "Error: could not open distances.csv\n";
        return 1;
    }

    // --- Process CSV with scaling ---
    // Parameters for scale_signal:
    // csv_path, column_index, csv_step, start_time, use_time,
    // observation_time, freq, integrator, check_drift, sec_per_burst
    
    int column_index = 0;
    double csv_step = (real_step > 0.0) ? real_step : 0.0001; // Default if not provided
    double start_time = 0.0;
    double use_time = simulation_time;
    double observation_time = simulation_time; // Use full time for observation
    double freq = 1.0 / csv_step; // Sampling frequency
    int integrator = RK4;
    bool check_drift = false; // Set to true if drift correction needed
    double sec_per_burst = -1.0; // Auto-detect
    
    ScaledSignalResult result = scale_signal(
        csv_file,
        column_index,
        csv_step,
        start_time,
        use_time,
        observation_time,
        freq,
        integrator,
        check_drift,
        sec_per_burst
    );
    
    if (!result.success) {
        std::cerr << "Error: scaling failed\n";
        return 1;
    }
    
    double dt_model = result.dt;

    // --- Initialize min/max and accumulated distance ---
    double h1_min = numeric_limits<double>::max();
    double h1_max = numeric_limits<double>::lowest();
    double syn_min = numeric_limits<double>::max();
    double syn_max = numeric_limits<double>::lowest();
    double accumulated_distance = 0.0;

    // --- Simulation loop ---
    for (size_t i = 0; i < result.scaled_signal.size(); ++i) {
        double time = i * dt_model;
        double csv_val = result.scaled_signal[i];

        // Overwrite dummy neuron
        HR_csv.set(HR::x, csv_val);

        double h1_val = h1.get(HR::x);
        s.step(dt_model, HR_csv.get(HR::x), h1_val);
        h1.add_synaptic_input(s.get(Synapsis::ifast));
        h1.step(dt_model);

        double syn_val = s.get(Synapsis::ifast);

        // Update min/max
        h1_min = min(h1_min, h1_val);
        h1_max = max(h1_max, h1_val);
        syn_min = min(syn_min, syn_val);
        syn_max = max(syn_max, syn_val);


        // Accumulate distance on-the-fly
        double h1_norm = (h1_val - h1_min) / (h1_max - h1_min + 1e-12);
        double syn_norm = (syn_val - syn_min) / (syn_max - syn_min + 1e-12);
        accumulated_distance += abs(h1_norm - syn_norm);

        // Write output
        out << time << " "
            << csv_val << " "
            << h1_val << " "
            << s.get(Synapsis::i) << " "
            << s.get(Synapsis::ifast) << " "
            << s.get(Synapsis::islow)
            << "\n";

    }
    std::cout << "H1  min: " << h1_min  << " | H1  max: " << h1_max  << std::endl;

    // --- Append distance to file ---
    if (!exists)
        dist_out << "timestamp,yaml_file,data_file,simulation_time,real_step,accumulated_distance\n";

    auto now = chrono::system_clock::now();
    time_t now_c = chrono::system_clock::to_time_t(now);

    dist_out << put_time(localtime(&now_c), "%F %T") << ","
             << yaml_file << ","
             << output_file << ","
             << simulation_time << ","
             << real_step << ","
             << accumulated_distance << "\n";

    out.close();
    dist_out.close();

    cout << "Simulation finished. Results written to " << output_file << "\n";
    cout << "Distance appended in " << dist_file << "\n";

    return 0;
}
