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

    std::string yaml_file   = argv[1];
    std::string output_file = argv[2];
    double simulation_time  = std::atof(argv[3]);
    double step = 0.0;
    std::string csv_file;
    if (argc == 6) {
        step = std::atof(argv[4]);
        csv_file = argv[5];
    } else {
        csv_file = argv[4];
    }

    cout << "YAML file: " << yaml_file << "\n";
    cout << "Output file: " << output_file << "\n";
    cout << "Simulation time: " << simulation_time << "\n";
    cout << "CSV file: " << csv_file << "\n";
    cout << "Step: " << step << "\n";

    // --- Initialize presynaptic neuron HR1 ---
    HR::ConstructorArgs args;
    auto param_names = HR::ParamNames();
    auto v_args = parse_file(yaml_file, "HR1");
    for(int i=0; i < HR::n_parameters; i++)
        args.params[i] = v_args[param_names[i]];

    HR h1(args);
    h1.set(HR::x, -0.712841);
    h1.set(HR::y, -1.93688);
    h1.set(HR::z, 3.16568);

    // --- Initialize dummy postsynaptic neuron for CSV ---
    HR::ConstructorArgs args_csv;
    HR HR_csv(args_csv);
    HR_csv.set(HR::x, 0.0); // initial value (will be overwritten by CSV)

    // --- Initialize synapse (HR1 -> CSV) ---
    Synapsis::ConstructorArgs syn_args;
    v_args = parse_file(yaml_file, "Chemical-HR1-HR2");
    auto syn_param_names = Synapsis::ParamNames();
    for(int i=0; i < Synapsis::n_parameters; i++)
        syn_args.params[i] = v_args[syn_param_names[i]];

    // HR1 presynaptic, HR_csv postsynaptic
    Synapsis s(h1, HR::x, HR_csv, HR::x, syn_args, 1);

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

    // --- Initialize min/max and accumulated distance ---
    double h1_min = numeric_limits<double>::max();
    double h1_max = numeric_limits<double>::lowest();
    double syn_min = numeric_limits<double>::max();
    double syn_max = numeric_limits<double>::lowest();
    double accumulated_distance = 0.0;
    double modelv_min = numeric_limits<double>::max();
    double modelv_max = numeric_limits<double>::lowest();
    double livingv_min = numeric_limits<double>::max();
    double livingv_max = numeric_limits<double>::lowest();

    size_t n_steps = static_cast<size_t>(simulation_time / step);

    // --- Simulation loop ---
    for (size_t i = 0; i < n_steps; ++i) {
        double time = i * step;
        double csv_val = 0.0;

        // Read next CSV line
        if (getline(csv, line)) {
            stringstream ss(line);
            string cell;
            if (getline(ss, cell, ',')) {
                try { csv_val = stod(cell); } catch (...) { csv_val = 0.0; }
            }
        }

        // Overwrite postsynaptic neuron with CSV value
        HR_csv.set(HR::x, csv_val);

        double h1_val = h1.get(HR::x)/1000;

        modelv_min = min(modelv_min, h1_val);
        modelv_max = max(modelv_max, h1_val);
        livingv_min = min(livingv_min, csv_val);
        livingv_max = max(livingv_max, csv_val);
        
        // Precalculated values to scale the HR to the recorded PD
        h1_val = 6.80447 * h1_val -0.0434337;

        // Step synapse
        s.step(step, h1_val, HR_csv.get(HR::x));

        // Apply synaptic input to postsynaptic neuron (optional, if you want to track input)
        HR_csv.add_synaptic_input(s.get(Synapsis::i));

        // Step presynaptic HR1 neuron
        h1.step(step);

        double syn_val = s.get(Synapsis::islow);

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
            << h1_val << " "    // presynaptic HR1
            << csv_val << " "   // postsynaptic CSV neuron
            << s.get(Synapsis::i) << " "
            << s.get(Synapsis::ifast) << " "
            << s.get(Synapsis::islow)
            << "\n";
    }

    csv.close();

    // --- Append distance to file ---
    if (!exists)
        dist_out << "timestamp,yaml_file,data_file,simulation_time,step,accumulated_distance\n";

    auto now = chrono::system_clock::now();
    time_t now_c = chrono::system_clock::to_time_t(now);

    dist_out << put_time(localtime(&now_c), "%F %T") << ","
             << yaml_file << ","
             << output_file << ","
             << simulation_time << ","
             << step << ","
             << accumulated_distance << "\n";

    out.close();
    dist_out.close();

    cout << "Simulation finished. Results written to " << output_file << "\n";
    cout << "Distance appended in " << dist_file << "\n";
    
    // Compute ranges
    double model_range = modelv_max - modelv_min;
    double living_range = livingv_max - livingv_min;
    
    // Compute factor and offset
    double factor_to_living = living_range / model_range;
    double offset_to_living = livingv_min - (modelv_min * factor_to_living);
    
    cout << "factor_to_living: " << factor_to_living << "\n";
    cout << "offset_to_living: " << offset_to_living << "\n";

    return 0;
}

