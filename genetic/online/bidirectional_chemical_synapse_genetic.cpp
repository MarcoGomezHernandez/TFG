/*
 * Copyright (C) 2011 Georgia Institute of Technology, University of Utah,
 * Weill Cornell Medical College
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "bidirectional_chemical_synapse_genetic.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <main_window.h>

extern "C" Plugin::Object *
createRTXIPlugin(void)
{
  return new BidirectionalChemicalSynapseGenetic();
}

static DefaultGUIModel::variable_t vars[] = {
    {"Individual evaluation time (s)", "Individual evaluation time; does not include stabilization time", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Individual stabilization time (s)", "Stabilization time for each individual; not included in Individual evaluation time", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Burst duration (s)", "-1 to use dynamic input", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Scale 2->1", "-1 to use dynamic input", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Offset 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Scale 1->2", "-1 to use dynamic input", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Offset 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Dynamic offset 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Dynamic offset 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},

    {"Neur 1 is living (1/0)", "Indicates if neuron 1 is alive (1) or Hindmarsh–Rose model (0)", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Neur 2 is living (1/0)", "Indicates if neuron 2 is alive (1) or Hindmarsh–Rose model (0)", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},

    // Config 2 -> 1
    {"E_syn 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_fast 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_fast 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_fast 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_slow 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k1 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k2 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_slow 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_slow 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Use I_fast 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_slow 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},

    // Config 1 -> 2
    {"E_syn 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_fast 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_fast 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_fast 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_slow 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k1 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k2 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_slow 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_slow 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Use I_fast 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_slow 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},

    {"Search phase genetic (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Current max to achieve genetic", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Current min to achieve genetic", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Current 2->1 (nA)", "Total synaptic current 2->1", DefaultGUIModel::OUTPUT},
    {"Current 1->2 (nA)", "Total synaptic current 1->2", DefaultGUIModel::OUTPUT},

    {"Voltage 1 (V)", "Membrane potential 1", DefaultGUIModel::INPUT},
    {"Voltage 2 (V)", "Membrane potential 2", DefaultGUIModel::INPUT},
    {"Scale 2->1", "Dynamic amplitude scale 2->1", DefaultGUIModel::INPUT},
    {"Offset 2->1", "Dynamic amplitude offset 2->1", DefaultGUIModel::INPUT},
    {"Scale 1->2", "Dynamic amplitude scale 1->2", DefaultGUIModel::INPUT},
    {"Offset 1->2", "Dynamic amplitude offset 1->2", DefaultGUIModel::INPUT},
    {"Burst duration 1 (s)", "Dynamic burst duration 1", DefaultGUIModel::INPUT},
    {"Burst duration 2 (s)", "Dynamic burst duration 2", DefaultGUIModel::INPUT},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

BidirectionalChemicalSynapseGenetic::BidirectionalChemicalSynapseGenetic(void)
    : DefaultGUIModel("RTHybrid Bidirectional Chemical Synapse Genetic", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>RTHybrid Bidirectional Chemical Synapse Genetic</b></p>");
  DefaultGUIModel::createGUI(vars, num_vars);
  customizeGUI();
  initParameters();
  update(INIT);
  refresh();
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

BidirectionalChemicalSynapseGenetic::~BidirectionalChemicalSynapseGenetic(void)
{
  if (genetic_NRT_thread.joinable())
  {
    stop_genetic.store(true, std::memory_order_relaxed);
    genetic_NRT_thread.join();
  }
}

void BidirectionalChemicalSynapseGenetic::runge_kutta_65(double (*f)(double, double, const ChemicalSynapseParams &), double &m_slow, double v_pre, double dt, const ChemicalSynapseParams &params)
{
  double apoyo, retorno;
  double k[6];

  retorno = (*f)(m_slow, v_pre, params);
  k[0] = dt * retorno;
  apoyo = m_slow + k[0] * 0.2;

  retorno = (*f)(apoyo, v_pre, params);
  k[1] = dt * retorno;
  apoyo = m_slow + k[0] * 0.075 + k[1] * 0.225;

  retorno = (*f)(apoyo, v_pre, params);
  k[2] = dt * retorno;
  apoyo = m_slow + k[0] * 0.3 - k[1] * 0.9 + k[2] * 1.2;

  retorno = (*f)(apoyo, v_pre, params);
  k[3] = dt * retorno;
  apoyo = m_slow + k[0] * 0.075 + k[1] * 0.675 - k[2] * 0.6 + k[3] * 0.75;

  retorno = (*f)(apoyo, v_pre, params);
  k[4] = dt * retorno;
  apoyo = m_slow + k[0] * 0.660493827160493 + k[1] * 2.5 - k[2] * 5.185185185185185 + k[3] * 3.888888888888889 - k[4] * 0.864197530864197;

  retorno = (*f)(apoyo, v_pre, params);
  k[5] = dt * retorno;

  m_slow += k[0] * 0.098765432098765 +
            k[2] * 0.396825396825396 +
            k[3] * 0.231481481481481 +
            k[4] * 0.308641975308641 -
            k[5] * 0.035714285714285;
}

void BidirectionalChemicalSynapseGenetic::select_dt_neuron_model(const double *dts, const double *pts, size_t length, double pts_live, double *dt, double *pts_burst)
{
  double aux = pts_live;
  double factor = 1;
  double intpart, fractpart;
  unsigned int flag = 0;
  int i;

  *dt = -1;
  *pts_burst = -1;

  while (aux < pts[0])
  {
    aux = pts_live * factor;
    factor += 1;

    for (i = length - 1; i >= 0; i--)
    {
      if (pts[i] > aux)
      {
        *dt = dts[i];
        *pts_burst = pts[i];
        fractpart = modf(*pts_burst / pts_live, &intpart);
        if (fractpart <= 0.1 * intpart)
          flag = 1;
        break;
      }
    }
    if (flag == 1)
      break;
  }

  if (flag == 0)
  {
    for (i = length - 1; i >= 0; i--)
    {
      if (pts[i] > aux)
      {
        *dt = dts[i];
        *pts_burst = pts[i];
        break;
      }
    }
  }
}

double BidirectionalChemicalSynapseGenetic::set_pts_burst(double sec_per_burst)
{
  constexpr size_t length = 144;
  const double pts_match = sec_per_burst * freq;
  double pts_burst;

  constexpr double dts[] = {0.000500, 0.000600, 0.000700, 0.000800, 0.000900, 0.001000, 0.001100, 0.001200, 0.001300, 0.001400, 0.001500, 0.001600, 0.001800, 0.002000, 0.002200, 0.002500, 0.002800, 0.002900, 0.003000, 0.003100, 0.003200, 0.003300, 0.003400, 0.003500, 0.003600, 0.003700, 0.003800, 0.003900, 0.004000, 0.004100, 0.004200, 0.004300, 0.004400, 0.004500, 0.004600, 0.004700, 0.004800, 0.004900, 0.005000, 0.005100, 0.005200, 0.005400, 0.005600, 0.005800, 0.006000, 0.006200, 0.006400, 0.006600, 0.006800, 0.007000, 0.007200, 0.007400, 0.007700, 0.008000, 0.008300, 0.008600, 0.008900, 0.009200, 0.009600, 0.010000, 0.010400, 0.010900, 0.011400, 0.011900, 0.012500, 0.013100, 0.013800, 0.014600, 0.015400, 0.016300, 0.017300, 0.018500, 0.019900, 0.021500, 0.023300, 0.025500, 0.028100, 0.028400, 0.028700, 0.029000, 0.029400, 0.029800, 0.030200, 0.030600, 0.031000, 0.031400, 0.031800, 0.032200, 0.032600, 0.033000, 0.033400, 0.033900, 0.034400, 0.034900, 0.035400, 0.035900, 0.036400, 0.036900, 0.037400, 0.038000, 0.038600, 0.039200, 0.039800, 0.040400, 0.041000, 0.041700, 0.042400, 0.043100, 0.043800, 0.044500, 0.045300, 0.046100, 0.046900, 0.047700, 0.048600, 0.049500, 0.050400, 0.051400, 0.052400, 0.053400, 0.054500, 0.055600, 0.056800, 0.058000, 0.059300, 0.060600, 0.062000, 0.063400, 0.064900, 0.066500, 0.068200, 0.069900, 0.071700, 0.073600, 0.075600, 0.077700, 0.079900, 0.082300, 0.084800, 0.087500, 0.090300, 0.093300, 0.096500, 0.100000};
  constexpr double pts[] = {577638.000000, 481366.000000, 412599.000000, 357615.500000, 317880.000000, 286092.500000, 259143.333333, 237548.000000, 218869.500000, 203236.000000, 189687.000000, 177634.000000, 157897.000000, 142001.833333, 129024.142857, 113496.125000, 101304.555556, 97811.222222, 94527.400000, 91478.200000, 88619.400000, 85916.636364, 83389.636364, 81007.090909, 78743.583333, 76615.416667, 74599.250000, 72676.000000, 70859.076923, 69130.846154, 67476.642857, 65907.357143, 64402.666667, 62971.466667, 61602.533333, 60286.187500, 59030.250000, 57825.562500, 56664.411765, 55553.294118, 54485.000000, 52463.222222, 50586.263158, 48841.842105, 47211.050000, 45685.666667, 44255.818182, 42914.772727, 41650.739130, 40459.083333, 39335.208333, 38270.680000, 36778.346154, 35398.000000, 34117.571429, 32926.517241, 31815.833333, 30777.612903, 29493.939394, 28313.588235, 27223.638889, 25974.405405, 24834.410256, 23790.268293, 22647.767442, 21609.977778, 20513.166667, 19388.627451, 18381.132075, 17365.719298, 16361.600000, 15299.937500, 14223.202899, 13164.400000, 12147.123457, 11098.876404, 10071.693878, 9965.282828, 9861.100000, 9759.059406, 9626.242718, 9497.009615, 9371.179245, 9248.672897, 9129.293578, 9012.981818, 8899.594595, 8789.000000, 8681.149123, 8575.896552, 8473.170940, 8348.176471, 8226.809917, 8108.934426, 7994.379032, 7883.015873, 7774.710938, 7669.348837, 7566.801527, 7447.308271, 7331.525926, 7219.291971, 7110.435714, 7004.816901, 6902.298611, 6786.417808, 6674.355705, 6565.940397, 6460.987013, 6359.339744, 6247.012579, 6138.592593, 6033.866667, 5932.660714, 5822.777778, 5716.896552, 5614.796610, 5505.541436, 5400.467391, 5299.324468, 5192.348958, 5089.615385, 4982.075000, 4878.985294, 4772.014354, 4669.633803, 4564.178899, 4463.381166, 4360.214912, 4255.294872, 4149.216667, 4048.296748, 3946.654762, 3844.764479, 3743.041353, 3641.872263, 3541.583630, 3438.300000, 3336.926421, 3233.951299, 3133.666667, 3032.899696, 2932.320588, 2829.684659};
  select_dt_neuron_model(dts, pts, length, pts_match, &dt, &pts_burst);

  return pts_burst;
}

double BidirectionalChemicalSynapseGenetic::sm_chemical_synapse_m(double m_slow, double v_pre, const ChemicalSynapseParams &params)
{
  return ((params.k1 * (1.0 - m_slow)) /
          (1.0 + exp(params.s_slow * (params.v_slow - v_pre)))) -
         (params.k2 * m_slow);
}

double BidirectionalChemicalSynapseGenetic::compute_i_slow(double &m_slow, double v_pre, double v_post, const ChemicalSynapseParams &params)
{
  for (size_t i = 0; i < s_points; i++)
  {
    runge_kutta_65(sm_chemical_synapse_m, m_slow, v_pre, dt, params);
  }
  return params.g_slow * m_slow * (v_post - params.e_syn);
}

double BidirectionalChemicalSynapseGenetic::compute_i_fast(double v_pre, double v_post, const ChemicalSynapseParams &params)
{
  return (params.g_fast * (v_post - params.e_syn)) /
         (1.0 + exp(params.s_fast * (params.v_fast - v_pre)));
}

void BidirectionalChemicalSynapseGenetic::execute(void)
{
  double val_i_slow_21 = 0.0, val_i_fast_21 = 0.0;
  double val_i_slow_12 = 0.0, val_i_fast_12 = 0.0;

  if (burst_duration_gui <= 0.0)
  {
    const double burst_duration_1 = input(6);
    const double burst_duration_2 = input(7);
    const double new_burst_duration = (burst_duration_1 == 0.0) ? burst_duration_2 : (burst_duration_2 == 0.0 ? burst_duration_1 : std::min(burst_duration_1, burst_duration_2));
    if (new_burst_duration != burst_duration)
    {
      burst_duration = new_burst_duration;
      s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq));
      if (s_points < 1)
        s_points = 1;
    }
  }

  const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;
  const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;

  double v1, v2;
  if (use_syn_21 || use_syn_12)
  {
    v1 = input(0) * 1000.0;
    v2 = input(1) * 1000.0;
  }

  size_t curr_scaling_factors_idx, new_scaling_factors_idx;
  double curr_scale, curr_offset;

  const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_acquire);

  double v2_scaled;
  if (use_syn_21)
  {
    curr_scaling_factors_idx = scaling_factors_21_idx.load(std::memory_order_relaxed);
    new_scaling_factors_idx = 1 - curr_scaling_factors_idx;
    curr_scale = scale_21[curr_scaling_factors_idx];
    curr_offset = offset_21[curr_scaling_factors_idx];
    double &new_scale = scale_21[new_scaling_factors_idx];
    double &new_offset = offset_21[new_scaling_factors_idx];

    new_scale = curr_scale;
    new_offset = curr_offset;

    if (scale_21_gui <= 0.0)
    {
      new_scale = input(2);
      if (new_scale == 0.0)
        new_scale = 1.0;
    }
    if (dynamic_offset_21)
      new_offset = input(3) * 1000.0;

    if (new_scale != curr_scale || new_offset != curr_offset)
      scaling_factors_21_idx.store(new_scaling_factors_idx, std::memory_order_release);

    v2_scaled = v2 * new_scale + new_offset;

    const ChemicalSynapseParams &curr_params_21 = params_21[curr_synapse_idx];
    if (use_i_slow_21)
      val_i_slow_21 = compute_i_slow(m_slow_21[curr_synapse_idx], v2_scaled, v1, curr_params_21);
    if (use_i_fast_21)
      val_i_fast_21 = compute_i_fast(v2_scaled, v1, curr_params_21);
  }

  double v1_scaled;
  if (use_syn_12)
  {
    curr_scaling_factors_idx = scaling_factors_12_idx.load(std::memory_order_relaxed);
    new_scaling_factors_idx = 1 - curr_scaling_factors_idx;
    curr_scale = scale_12[curr_scaling_factors_idx];
    curr_offset = offset_12[curr_scaling_factors_idx];
    double &new_scale = scale_12[new_scaling_factors_idx];
    double &new_offset = offset_12[new_scaling_factors_idx];

    new_scale = curr_scale;
    new_offset = curr_offset;

    if (scale_12_gui <= 0.0)
    {
      new_scale = input(2);
      if (new_scale == 0.0)
        new_scale = 1.0;
    }
    if (dynamic_offset_12)
      new_offset = input(3) * 1000.0;

    if (new_scale != curr_scale || new_offset != curr_offset)
      scaling_factors_12_idx.store(new_scaling_factors_idx, std::memory_order_release);

    v1_scaled = v1 * new_scale + new_offset;

    const ChemicalSynapseParams &curr_params_12 = params_12[curr_synapse_idx];
    if (use_i_slow_12)
      val_i_slow_12 = compute_i_slow(m_slow_12[curr_synapse_idx], v1_scaled, v2, curr_params_12);
    if (use_i_fast_12)
      val_i_fast_12 = compute_i_fast(v1_scaled, v2, curr_params_12);
  }

  if (RT_storing.load(std::memory_order_acquire))
  {
    if (storing_idx < num_elements)
    {
      if (use_syn_21)
      {
        v1_sig[storing_idx] = v1;
        v2_scaled_sig[storing_idx] = v2_scaled;
        if (use_i_fast_21)
          i_fast_sig_21[storing_idx] = val_i_fast_21;
        if (use_i_slow_21)
          i_slow_21[storing_idx] = val_i_slow_21;
      }
      if (use_syn_12)
      {
        v2_sig[storing_idx] = v2;
        v1_scaled_sig[storing_idx] = v1_scaled;
        if (use_i_fast_12)
          i_fast_sig_12[storing_idx] = val_i_fast_12;
        if (use_i_slow_12)
          i_slow_12[storing_idx] = val_i_slow_12;
      }
      storing_idx++;
    }
    else
    {
      RT_storing.store(false, std::memory_order_release);
    }
  }

  output(0) = val_i_fast_21 + val_i_slow_21;
  output(1) = val_i_fast_12 + val_i_slow_12;
}

void BidirectionalChemicalSynapseGenetic::initParameters(void)
{
  evaluation_time = 10.0;
  stabilization_time = 1.0;

  search_phase = 1u;
  current_max = 10.0;
  current_min = -10.0;

  burst_duration_gui = 1.0;
  burst_duration = burst_duration_gui;

  is_living_1 = 0u;
  is_living_2 = 0u;
  dynamic_offset_21 = 0u;
  dynamic_offset_12 = 0u;
  scale_21_gui = 1.0;
  scale_21[0] = scale_21_gui;
  offset_21[0] = 0.0 * 1000.0;
  scale_12_gui = 1.0;
  scale_12[0] = scale_12_gui;
  offset_12[0] = 0.0 * 1000.0;

  m_slow_21[0] = 0.0;
  m_slow_12[0] = 0.0;

  stop_genetic.store(false, std::memory_order_relaxed);
  RT_storing.store(false, std::memory_order_relaxed);
  genetic_running.store(false, std::memory_order_relaxed);

  synapse_idx.store(0, std::memory_order_relaxed);
  scaling_factors_21_idx.store(0, std::memory_order_relaxed);
  scaling_factors_12_idx.store(0, std::memory_order_relaxed);

  params_21[0].e_syn = -1.92;
  params_21[0].g_fast = 0.046;
  params_21[0].s_fast = 0.44;
  params_21[0].v_fast = -1.66;
  params_21[0].g_slow = 0.208;
  params_21[0].k1 = 0.74;
  params_21[0].k2 = 0.007;
  params_21[0].s_slow = 1.0;
  params_21[0].v_slow = -1.74;

  use_i_fast_21 = 1u;
  use_i_slow_21 = 1u;

  params_12[0].e_syn = -1.92;
  params_12[0].g_fast = 0.046;
  params_12[0].s_fast = 0.44;
  params_12[0].v_fast = -1.66;
  params_12[0].g_slow = 0.208;
  params_12[0].k1 = 0.74;
  params_12[0].k2 = 0.007;
  params_12[0].s_slow = 1.0;
  params_12[0].v_slow = -1.74;

  use_i_fast_12 = 1u;
  use_i_slow_12 = 1u;
}

void BidirectionalChemicalSynapseGenetic::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag)
  {
  case INIT:
    period = RT::System::getInstance()->getPeriod() * 1e-9; // s
    freq = 1.0 / period;
    s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq));
    if (s_points == 0)
      s_points = 1;

    setParameter("Individual evaluation time (s)", evaluation_time);
    setParameter("Individual stabilization time (s)", stabilization_time);
    setParameter("Burst duration (s)", burst_duration_gui);

    setParameter("Search phase genetic (1/0)", search_phase);
    setParameter("Current max to achieve genetic", current_max);
    setParameter("Current min to achieve genetic", current_min);

    setParameter("Neur 1 is living (1/0)", is_living_1);
    setParameter("Neur 2 is living (1/0)", is_living_2);
    setParameter("Scale 2->1", scale_21_gui);
    setParameter("Offset 2->1", offset_21[scaling_factors_21_idx.load(std::memory_order_relaxed)] / 1000.0);
    setParameter("Scale 1->2", scale_12_gui);
    setParameter("Offset 1->2", offset_12[scaling_factors_12_idx.load(std::memory_order_relaxed)] / 1000.0);

    setParameter("Dynamic offset 2->1 (1/0)", dynamic_offset_21);
    setParameter("Dynamic offset 1->2 (1/0)", dynamic_offset_12);

    update_params_gui();

    setParameter("Use I_fast 2->1 (1/0)", use_i_fast_21);
    setParameter("Use I_slow 2->1 (1/0)", use_i_slow_21);
    setParameter("Use I_fast 1->2 (1/0)", use_i_fast_12);
    setParameter("Use I_slow 1->2 (1/0)", use_i_slow_12);

    break;

  case MODIFY:
  {
    burst_duration_gui = getParameter("Burst duration (s)").toDouble();
    if ((burst_duration_gui > 0.0) && (burst_duration_gui != burst_duration))
    {
      burst_duration = burst_duration_gui;
      s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq));
      if (s_points < 1)
        s_points = 1;
    }

    size_t curr_scaling_factors_idx, new_scaling_factors_idx;
    double curr_scale, curr_offset;

    scale_21_gui = getParameter("Scale 2->1").toDouble();
    dynamic_offset_21 = getParameter("Dynamic offset 2->1 (1/0)").toUInt();
    curr_scaling_factors_idx = scaling_factors_21_idx.load(std::memory_order_relaxed);
    new_scaling_factors_idx = 1 - curr_scaling_factors_idx;
    curr_scale = scale_21[curr_scaling_factors_idx];
    curr_offset = offset_21[curr_scaling_factors_idx];
    {
      double &new_scale = scale_21[new_scaling_factors_idx];
      double &new_offset = offset_21[new_scaling_factors_idx];

      new_scale = curr_scale;
      new_offset = curr_offset;

      if (scale_21_gui > 0.0)
        new_scale = scale_21_gui;
      if (!dynamic_offset_21)
        new_offset = getParameter("Offset 2->1").toDouble() * 1000.0;

      if (new_scale != curr_scale || new_offset != curr_offset)
        scaling_factors_21_idx.store(new_scaling_factors_idx, std::memory_order_release);
    }

    scale_12_gui = getParameter("Scale 1->2").toDouble();
    dynamic_offset_12 = getParameter("Dynamic offset 1->2 (1/0)").toUInt();
    curr_scaling_factors_idx = scaling_factors_12_idx.load(std::memory_order_relaxed);
    new_scaling_factors_idx = 1 - curr_scaling_factors_idx;
    curr_scale = scale_12[curr_scaling_factors_idx];
    curr_offset = offset_12[curr_scaling_factors_idx];
    {
      double &new_scale = scale_12[new_scaling_factors_idx];
      double &new_offset = offset_12[new_scaling_factors_idx];

      new_scale = curr_scale;
      new_offset = curr_offset;

      if (scale_12_gui > 0.0)
        new_scale = scale_12_gui;
      if (!dynamic_offset_12)
        new_offset = getParameter("Offset 1->2").toDouble() * 1000.0;

      if (new_scale != curr_scale || new_offset != curr_offset)
        scaling_factors_12_idx.store(new_scaling_factors_idx, std::memory_order_release);
    }

    if (!genetic_running.load(std::memory_order_acquire))
    {
      evaluation_time = getParameter("Individual evaluation time (s)").toDouble();
      stabilization_time = getParameter("Individual stabilization time (s)").toDouble();
      is_living_1 = getParameter("Neur 1 is living (1/0)").toUInt();
      is_living_2 = getParameter("Neur 2 is living (1/0)").toUInt();

      search_phase = getParameter("Search phase genetic (1/0)").toUInt();
      current_max = getParameter("Current max to achieve genetic").toDouble();
      current_min = getParameter("Current min to achieve genetic").toDouble();

      const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_acquire);

      params_21[curr_synapse_idx].e_syn = getParameter("E_syn 2->1").toDouble();
      params_21[curr_synapse_idx].g_fast = getParameter("g_fast 2->1").toDouble();
      params_21[curr_synapse_idx].s_fast = getParameter("s_fast 2->1").toDouble();
      params_21[curr_synapse_idx].v_fast = getParameter("V_fast 2->1").toDouble();
      params_21[curr_synapse_idx].g_slow = getParameter("g_slow 2->1").toDouble();
      params_21[curr_synapse_idx].k1 = getParameter("k1 2->1").toDouble();
      params_21[curr_synapse_idx].k2 = getParameter("k2 2->1").toDouble();
      params_21[curr_synapse_idx].s_slow = getParameter("s_slow 2->1").toDouble();
      params_21[curr_synapse_idx].v_slow = getParameter("V_slow 2->1").toDouble();
      use_i_fast_21 = getParameter("Use I_fast 2->1 (1/0)").toUInt();
      use_i_slow_21 = getParameter("Use I_slow 2->1 (1/0)").toUInt();

      params_12[curr_synapse_idx].e_syn = getParameter("E_syn 1->2").toDouble();
      params_12[curr_synapse_idx].g_fast = getParameter("g_fast 1->2").toDouble();
      params_12[curr_synapse_idx].s_fast = getParameter("s_fast 1->2").toDouble();
      params_12[curr_synapse_idx].v_fast = getParameter("V_fast 1->2").toDouble();
      params_12[curr_synapse_idx].g_slow = getParameter("g_slow 1->2").toDouble();
      params_12[curr_synapse_idx].k1 = getParameter("k1 1->2").toDouble();
      params_12[curr_synapse_idx].k2 = getParameter("k2 1->2").toDouble();
      params_12[curr_synapse_idx].s_slow = getParameter("s_slow 1->2").toDouble();
      params_12[curr_synapse_idx].v_slow = getParameter("V_slow 1->2").toDouble();
      use_i_fast_12 = getParameter("Use I_fast 1->2 (1/0)").toUInt();
      use_i_slow_12 = getParameter("Use I_slow 1->2 (1/0)").toUInt();
    }

    break;
  }

  case UNPAUSE:
    break;

  case PERIOD:
  {
    const double new_period = RT::System::getInstance()->getPeriod() * 1e-9; // s
    if (new_period != period)
    {
      period = new_period;
      freq = 1.0 / period;
      s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq));
      if (s_points == 0)
        s_points = 1;

      if (genetic_running.load(std::memory_order_relaxed))
      {
        stop_genetic.store(true, std::memory_order_relaxed); // Porque cambiaría el número de puntos a almacenar
      }
    }

    break;
  }

  case PAUSE:
    output(0) = 0;
    output(1) = 0;
    break;

  default:
    break;
  }
}

void BidirectionalChemicalSynapseGenetic::customizeGUI(void)
{
  QGridLayout *customlayout = DefaultGUIModel::getLayout();
  gentic_button = new QPushButton("Start Genetic");
  QObject::connect(gentic_button, SIGNAL(clicked()), this, SLOT(toggle_genetic_event()));
  customlayout->addWidget(gentic_button, 0, 0);
  setLayout(customlayout);
}

void BidirectionalChemicalSynapseGenetic::toggle_genetic_event(void)
{
  if (!genetic_running.load(std::memory_order_relaxed))
  {
    // Lógica para EMPEZAR
    if (genetic_NRT_thread.joinable())
    {
      genetic_NRT_thread.join();
      stop_genetic.store(false, std::memory_order_relaxed);
    }
    genetic_running.store(true, std::memory_order_relaxed);
    gentic_button->setText("Stop Genetic");
    set_params_read_only(true);
    size_t curr_scaling_factors_21_idx = scaling_factors_21_idx.load(std::memory_order_relaxed);
    size_t curr_scaling_factors_12_idx = scaling_factors_12_idx.load(std::memory_order_relaxed);
    genetic_NRT_thread = std::thread(&BidirectionalChemicalSynapseGenetic::NRT_genetic, this, freq, scale_21[curr_scaling_factors_21_idx], offset_21[curr_scaling_factors_21_idx], scale_12[curr_scaling_factors_12_idx], offset_12[curr_scaling_factors_12_idx]);
  }
  else
  {
    // Lógica para PARAR
    stop_genetic.store(true, std::memory_order_relaxed);
  }
}

void BidirectionalChemicalSynapseGenetic::stop_genetic_event_async(void)
{
  set_params_read_only(false);
  gentic_button->setText("Start Genetic");
  genetic_running.store(false, std::memory_order_release);
}

void BidirectionalChemicalSynapseGenetic::set_params_read_only(bool read_only)
{
  DefaultGUILineEdit *example_edit = parameter["Individual evaluation time (s)"].edit;
  QPalette palette = example_edit->palette;
  palette.setBrush(example_edit->foregroundRole(), read_only ? Qt::darkGray : QApplication::palette().color(QPalette::WindowText));

  set_param_read_only("Individual evaluation time (s)", palette, read_only);
  set_param_read_only("Individual stabilization time (s)", palette, read_only);
  set_param_read_only("Neur 1 is living (1/0)", palette, read_only);
  set_param_read_only("Neur 2 is living (1/0)", palette, read_only);

  set_param_read_only("Search phase genetic (1/0)", palette, read_only);
  set_param_read_only("Current max to achieve genetic", palette, read_only);
  set_param_read_only("Current min to achieve genetic", palette, read_only);

  set_param_read_only("E_syn 2->1", palette, read_only);
  set_param_read_only("g_fast 2->1", palette, read_only);
  set_param_read_only("s_fast 2->1", palette, read_only);
  set_param_read_only("V_fast 2->1", palette, read_only);
  set_param_read_only("g_slow 2->1", palette, read_only);
  set_param_read_only("k1 2->1", palette, read_only);
  set_param_read_only("k2 2->1", palette, read_only);
  set_param_read_only("s_slow 2->1", palette, read_only);
  set_param_read_only("V_slow 2->1", palette, read_only);
  set_param_read_only("Use I_fast 2->1 (1/0)", palette, read_only);
  set_param_read_only("Use I_slow 2->1 (1/0)", palette, read_only);

  set_param_read_only("E_syn 1->2", palette, read_only);
  set_param_read_only("g_fast 1->2", palette, read_only);
  set_param_read_only("s_fast 1->2", palette, read_only);
  set_param_read_only("V_fast 1->2", palette, read_only);
  set_param_read_only("g_slow 1->2", palette, read_only);
  set_param_read_only("k1 1->2", palette, read_only);
  set_param_read_only("k2 1->2", palette, read_only);
  set_param_read_only("s_slow 1->2", palette, read_only);
  set_param_read_only("V_slow 1->2", palette, read_only);
  set_param_read_only("Use I_fast 1->2 (1/0)", palette, read_only);
  set_param_read_only("Use I_slow 1->2 (1/0)", palette, read_only);
}

void BidirectionalChemicalSynapseGenetic::set_param_read_only(const QString &name, const QPalette &pal, bool read_only)
{
  DefaultGUILineEdit *edit = parameter[name].edit;
  edit->setReadOnly(read_only);
  edit->setPalette(pal);
}

void BidirectionalChemicalSynapseGenetic::update_params_gui(void)
{
  const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);

  setParameter("E_syn 2->1", params_21[curr_synapse_idx].e_syn);
  setParameter("g_fast 2->1", params_21[curr_synapse_idx].g_fast);
  setParameter("s_fast 2->1", params_21[curr_synapse_idx].s_fast);
  setParameter("V_fast 2->1", params_21[curr_synapse_idx].v_fast);
  setParameter("g_slow 2->1", params_21[curr_synapse_idx].g_slow);
  setParameter("k1 2->1", params_21[curr_synapse_idx].k1);
  setParameter("k2 2->1", params_21[curr_synapse_idx].k2);
  setParameter("s_slow 2->1", params_21[curr_synapse_idx].s_slow);
  setParameter("V_slow 2->1", params_21[curr_synapse_idx].v_slow);

  setParameter("E_syn 1->2", params_12[curr_synapse_idx].e_syn);
  setParameter("g_fast 1->2", params_12[curr_synapse_idx].g_fast);
  setParameter("s_fast 1->2", params_12[curr_synapse_idx].s_fast);
  setParameter("V_fast 1->2", params_12[curr_synapse_idx].v_fast);
  setParameter("g_slow 1->2", params_12[curr_synapse_idx].g_slow);
  setParameter("k1 1->2", params_12[curr_synapse_idx].k1);
  setParameter("k2 1->2", params_12[curr_synapse_idx].k2);
  setParameter("s_slow 1->2", params_12[curr_synapse_idx].s_slow);
  setParameter("V_slow 1->2", params_12[curr_synapse_idx].v_slow);
}