#ifndef FITNESS_ONLINE_H
#define FITNESS_ONLINE_H

#include <array>
#include <kfr/all.hpp>
using namespace kfr;

struct ChemicalSynapseVariationParams
{
    double g_fast;
    double g_slow;
    double s_fast;
    double v_fast;
    double v_slow;
    double k1;
    double k2;
    double s_slow;
};

struct Individual
{
  ChemicalSynapseVariationParams params;
  double fitness;
};

struct FitnessExtraData
{
  univector<double> ifast_sig_1;
  univector<double> islow_sig_1;
  univector<double> ifast_sig_2;
  univector<double> islow_sig_2;

  bool phase;
  bool antiphase;

  double fs;
  double filter_fc;

  double i_range_expected_min_syn1;
  double i_range_expected_max_syn1;
  double i_range_expected_min_syn2;
  double i_range_expected_max_syn2;

  bool use_ifast;
  bool use_islow;
};

struct SigBuffers
{

  univector<double> v1_sig;
  univector<double> v2_sig;
  univector<double> syn1_ifast_sig;
  univector<double> syn1_islow_sig;
  univector<double> syn2_ifast_sig;
  univector<double> syn2_islow_sig;


  univector<double> ref_ifast_centered;
  univector<double> ref_islow_centered;
};

double fitness_from_sigs(const FitnessExtraData &extra,
                         SigBuffers &buffers);

#endif
