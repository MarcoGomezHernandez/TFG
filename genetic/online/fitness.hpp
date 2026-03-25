#ifndef FITNESS_H
#define FITNESS_H

#include <array>
#include <kfr/all.hpp>
using namespace kfr;

struct ChemicalSynapsisVariationParams
{
    double gfast;
    double gslow;
    double sfast;
    double Vfast;
    double Vslow;
    double k1;
    double k2;
    double sslow;
};

struct Individual
{
  ChemicalSynapsisVariationParams params;
  double fitness;
};

struct ConstantSigFitnessVals
{
  univector<double> ifast_sig_centered_to_fit;
  double ifast_sig_factor_to_fit;
  univector<double> islow_sig_centered_to_fit;
  double islow_sig_factor_to_fit;

  double ifast_sig_min_to_fit;
  double ifast_sig_max_to_fit;
  double ifast_sig_range_to_fit;

  double islow_sig_min_to_fit;
  double islow_sig_max_to_fit;
  double islow_sig_range_to_fit;
};

struct SigBuffers
{
  univector<double> ifast_sig;
  univector<double> islow_sig;
};

ConstantSigFitnessVals calc_const_sig_fitness_vals(const univector_ref<const double> &vpre_sig,
                                                   double pts_burst_real,
                                                   bool use_ifast,
                                                   bool use_islow,
                                                   bool search_phase);

#endif
