/*************************************************************

Copyright (c) 2006, Fernando Herrero Carrón
Copyright (c) 2020, Angel Lareo <angel.lareo@gmail.com>
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

#include <DifferentialNeuronWrapper.h>
#include <HindmarshRoseModel.h>
#include <RungeKutta4.h>
#include <SystemWrapper.h>
#include <iostream>

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<HindmarshRoseModel<double>>, Integrator>
    Neuron;

int main(int argc, char **argv) {
  // Struct to initialize neuron model parameters
  Neuron::ConstructorArgs args;

  // Set the parameter values
  args.params[Neuron::e] = 3.281;
  args.params[Neuron::mu] = 0.0029;
  args.params[Neuron::S] = 4;
  args.params[Neuron::a] = 1;
  args.params[Neuron::b] = 3;
  args.params[Neuron::c] = 1;
  args.params[Neuron::d] = 5;
  args.params[Neuron::xr] = -1.6;
  args.params[Neuron::vh] = 1;

  // Initialize a new neuron model
  Neuron n(args);

  // You can also initialize the variables of the neuron model to a given value
  n.set(Neuron::x, -0.712841);
  n.set(Neuron::y, -1.93688);
  n.set(Neuron::z, 3.16568);

  // Set the integration step
  const double step = 0.01;

  // Perform the simulation
  double simulation_time = 10000;

  std::cout << "time " << "HR1" << std::endl;

  int flag=0;
  for (double time = 0; time < simulation_time; time += step) {
    n.step(step);

    if (flag%3)
      std::cout << time << " " << n.get(Neuron::x) << std::endl;    

    flag++;
  }

  return 0;
}
