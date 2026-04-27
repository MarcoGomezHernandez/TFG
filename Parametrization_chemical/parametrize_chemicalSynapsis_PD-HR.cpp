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

    double model_step = 0.01;
    int step_factor = 8;

    cout << "YAML file: " << yaml_file << "\n";
    cout << "Output file: " << output_file << "\n";
    cout << "Simulation time: " << simulation_time << "\n";
    cout << "CSV file: " << csv_file << "\n";
    cout << "Real step: " << real_step << "\n";
    cout << "Model step: " << model_step << "\n";
    cout << "Step factor: " << step_factor << "\n";

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

    // --- Open CSV file for on-the-fly reading ---
    ifstream csv(csv_file);
    if (!csv.is_open()) throw runtime_error("Could not open csv_file");
    string line;

    // --- Pre-pass to compute livingv_min and livingv_max ---
    double livingv_min = numeric_limits<double>::max();
    double livingv_max = numeric_limits<double>::lowest();
    while (getline(csv, line)) {
        stringstream ss(line);
        string cell;
        if (getline(ss, cell, ',')) {
            try {
                double val = stod(cell);
                livingv_min = min(livingv_min, val);
                livingv_max = max(livingv_max, val);
            } catch (...) {
                // Skip invalid lines
            }
        }
    }
    // Reset file pointer to beginning
    csv.clear();
    csv.seekg(0);

    // --- Compute scaling factors ---
    double modelv_min = -1.608734;
    double modelv_max = 1.797032;
    double model_range = modelv_max - modelv_min;
    double living_range = livingv_max - livingv_min;
    double factor_to_model = model_range / living_range;
    double offset_to_model = modelv_min - (livingv_min * factor_to_model);

    // --- Initialize min/max and accumulated distance ---
    double h1_min = numeric_limits<double>::max();
    double h1_max = numeric_limits<double>::lowest();
    double syn_min = numeric_limits<double>::max();
    double syn_max = numeric_limits<double>::lowest();
    double accumulated_distance = 0.0;

    size_t n_steps = static_cast<size_t>(simulation_time / real_step);
    

    // --- Simulation loop ---
    for (size_t i = 0; i < n_steps; ++i) {
        double time = i * real_step;
        double csv_val = 0.0;

        // Read next CSV line
        if (getline(csv, line)) {
            stringstream ss(line);
            string cell;
            if (getline(ss, cell, ',')) {
                try { csv_val = stod(cell); } catch (...) { csv_val = 0.0; }
            }
        }

        // Apply scaling to csv_val
        csv_val = factor_to_model * csv_val + offset_to_model;

        // Overwrite dummy neuron
        HR_csv.set(HR::x, csv_val);

        double h1_val = 0.0;
        for (int j = 0; j < step_factor; ++j) {
            h1_val = h1.get(HR::x);
            s.step(model_step, HR_csv.get(HR::x), h1_val);
            h1.add_synaptic_input(-s.get(Synapsis::i));
            h1.step(model_step);
        }

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
    csv.close();

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
