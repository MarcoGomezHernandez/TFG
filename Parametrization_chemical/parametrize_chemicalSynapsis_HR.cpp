/*************************************************************

Copyright (c) 2025, Alicia Garrido Peña <alicia.garrido@uam.es>
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
#include <chrono>
#include <iomanip>
#include <ctime>
#include <filesystem>

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
// typedef ChemicalSynapsisModel<double> SynapsisModel;


std::map<std::string,double> parse_file(const std::string& file, const std::string& section)
{
    std::map<std::string,double> config_values;

    // Cargar el YAML
    YAML::Node config = YAML::LoadFile(file);

    // Verificar que exista la subsección
    if (!config[section]) {
        std::cerr << "Section '" << section << "' not found in " << file << std::endl;
        return config_values;
    }

    // Iterar por los pares clave/valor de la subsección
    for (auto it = config[section].begin(); it != config[section].end(); ++it) {
        std::string key = it->first.as<std::string>();
        double value = it->second.as<double>();
        config_values[key] = value;
    }

    return config_values;
}


int main(int argc, char **argv) {
  if (argc < 4) {
      std::cerr << "Usage: " << argv[0]
                << " <yaml_file> <output_file> <simulation_time> [step]\n";
      return 1;
  }

  // Required arguments
  std::string yaml_file   = argv[1];
  std::string output_file = argv[2];
  double simulation_time  = std::atof(argv[3]);

  // Optional argument (step)
  double step = 0.0; // default value
  if (argc >= 5) {
      step = std::atof(argv[4]);
  }

  // Echo input back
  std::cout << "YAML file: " << yaml_file << "\n";
  std::cout << "Output file: " << output_file << "\n";
  std::cout << "Simulation time: " << simulation_time << "\n";
  std::cout << "Step: " << step << "\n";

  // Initialize neuron model parameters
  HR::ConstructorArgs args, args2; 

  std::vector<std::string> param_names = HR::ParamNames();

  std::map<std::string,double> v_args = parse_file(yaml_file, "HR1");

  for(int i=0; i < HR::n_parameters; i++)
  {
    // std::cout << v_args[param_names[i]] << " " << param_names[i] << endl;
    args.params[i] = v_args[param_names[i]];
  }

  HR h1(args);
  h1.set(HR::x, -0.712841);
  h1.set(HR::y, -1.93688);
  h1.set(HR::z, 3.16568);

  v_args = parse_file(yaml_file, "HR2");

  for(int i=0; i < HR::n_parameters; i++)
  {
    // std::cout << v_args[param_names[i]] << " " << param_names[i] << endl;
    args2.params[i] = v_args[param_names[i]];
  }

  HR h2(args2);
  
  h2.set(HR::x, -0.712841);
  h2.set(HR::y, -1.93688);
  h2.set(HR::z, 3.16568);

  Synapsis::ConstructorArgs syn_args;
  v_args = parse_file(yaml_file, "Chemical-HR1-HR2");
  param_names = Synapsis::ParamNames();

  for(int i=0; i < Synapsis::n_parameters; i++)
  {
    // std::cout << v_args[param_names[i]] << " " << param_names[i] << endl;
    syn_args.params[i] = v_args[param_names[i]];
  }
  // Initialize a synapsis between the neurons
  Synapsis s(h1, HR::x, h2, HR::x, syn_args, 1);
  

  // Set the parameter values
  /*args.params[Neuron::e] = 3.281;
  args.params[Neuron::mu] = 0.0029;
  args.params[Neuron::S] = 4;
  args.params[Neuron::a] = 1;
  args.params[Neuron::b] = 3;
  args.params[Neuron::c] = 1;
  args.params[Neuron::d] = 5;
  args.params[Neuron::xr] = -1.6;
  args.params[Neuron::vh] = 1;

  Synapsis::ConstructorArgs syn_args;
  syn_args.params[Synapsis::gfast] = 0.015;
  syn_args.params[Synapsis::Esyn] = -75;
  syn_args.params[Synapsis::sfast] = 0.2;
  syn_args.params[Synapsis::Vfast] = -50;
  syn_args.params[Synapsis::gslow] = 0.025; //When 0, use only fast
  syn_args.params[Synapsis::k1] = 1;
  syn_args.params[Synapsis::k2] = 0.03;
  syn_args.params[Synapsis::sslow] = 1;*/


  // // Set initial value of V in neuron n1
  // h1.set(HR::v, -75);

  // Open output file
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Error: could not open " << output_file << " for writing.\n";
        return 1;
    }

    // Write header
    out << "Time Vpre Vpost i ifast islow\n";

    // Check if distances.txt exists already
    std::string dist_file = "distances.csv";
    bool exists = std::filesystem::exists(dist_file);

    // Secondary file (append mode)
    std::ofstream dist_out(dist_file, std::ios::app);
    if (!dist_out.is_open()) {
        std::cerr << "Error: could not open distances.csv\n";
        return 1;
    }

    double h1_min = std::numeric_limits<double>::max();
    double h1_max = std::numeric_limits<double>::lowest();
    double syn_min = std::numeric_limits<double>::max();
    double syn_max = std::numeric_limits<double>::lowest();

    // Option 1: store all values in vectors for second pass
    std::vector<double> h1_vals;
    std::vector<double> syn_vals;
    std::vector<double> times;

    // Simulation loop
    for (double time = 0; time < simulation_time; time += step) {
        s.step(step, h1.get(HR::x), h2.get(HR::x));

        // Provide an external current input to both neurons
        h2.add_synaptic_input(s.get(Synapsis::i));

        h1.step(step);
        h2.step(step);

        double h1_val = h1.get(HR::x);
        double syn_val = s.get(Synapsis::islow);

        // store
        times.push_back(time);
        h1_vals.push_back(h1_val);
        syn_vals.push_back(syn_val);

        // update min/max
        if (h1_val < h1_min) h1_min = h1_val;
        if (h1_val > h1_max) h1_max = h1_val;
        if (syn_val < syn_min) syn_min = syn_val;
        if (syn_val > syn_max) syn_max = syn_val;

        out << time << " "
            << h1.get(HR::x) << " "
            << h2.get(HR::x) << " "
            << s.get(Synapsis::i) << " "
            << s.get(Synapsis::ifast) << " "
            << s.get(Synapsis::islow)
            << "\n";
    }

    double accumulated_distance = 0.0;
    for (size_t i = 0; i < h1_vals.size(); ++i) {
        double h1_norm = (h1_vals[i] - h1_min) / (h1_max - h1_min);
        double syn_norm = (syn_vals[i] - syn_min) / (syn_max - syn_min);

        double distance = std::abs(h1_norm - syn_norm);
        accumulated_distance += distance;
    } 

    // Write header if first time
    if (!exists) {
        dist_out << "timestamp,yaml_file,data_file,simulation_time,step,accumulated_distance\n";
    }
    // Get current time
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    // After loop: append just once
    dist_out << std::put_time(std::localtime(&now_c), "%F %T") << "," 
            << yaml_file << ","
            << output_file << ","
            << simulation_time << ","
            << step << ","
            << accumulated_distance << "\n";

    out.close();
    dist_out.close();
    std::cout << "Simulation finished. Results written to " << output_file << "\n";
    std::cout << "Distance append in " << dist_file << "\n";


  return 0;
}
