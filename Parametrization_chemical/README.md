This is a repository for an automated parametrization of the chemical synapse by Gollowasch et al., between a Hindmarsh-Rose model neuron and a recorded PD from the pyloric CPG (crustacean STG).

Installation:
1. Install NeuN library from https://github.com/GNB-UAM/neun
2. Clone this repository

Usage:
1. yaml_generator.py: It will generate all possible combinations of parameters in yaml_configs
2. .cpp files: would run the specified circuit of model/recorded neurons and synapses between them
3. run_parametrization.py: It will run the parametrization of the desired setup specified in a .cpp with all the configs inside yaml_configs

Inside the .cpp it can be specified to save a desired value into a file to keep track of how each configuration has performed to be analyzed automatically. For example saving the distance between the neurons to know how effective each configuration of parameters has been in synchronyzying the neurons
