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
    // Genetic parameters
    {"Genetic generations completed", "", DefaultGUIModel::STATE},
    {"Genetic individuals of the generation completed", "", DefaultGUIModel::STATE},
    {"Genetic num generations", "Number of generations for the genetic algorithm", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Genetic population size", "Population size for the genetic algorithm", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Individual evaluation time genetic (s)", "Individual evaluation time; does not include stabilization time", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Individual stabilization time genetic (s)", "Stabilization time for each individual; not included in Individual evaluation time", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Search phase genetic (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Current ranges to achieve are from scale of neuron (1/2)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Current max to achieve genetic 2->1", "Both directions are in the same scale", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Current min to achieve genetic 2->1", "Both directions are in the same scale", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Current max to achieve genetic 1->2", "Both directions are in the same scale", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Current min to achieve genetic 1->2", "Both directions are in the same scale", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    // Use fast/slow currents
    {"Use I_fast 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_slow 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_fast 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_slow 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},

    // Dynamic, Offsets, Scales, Min/Max
    {"Dynamic offsets (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Dynamic min and max 1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Offset 1->2", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Offset 2->1", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Scale 1->2", "-1 to use dynamic input", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Scale 2->1", "-1 to use dynamic input", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Max 1 (V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Min 1 (V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"dt", "Integration timestep for synapse model", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

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

    {"Current 2->1 (nA)", "Total synaptic current 2->1", DefaultGUIModel::OUTPUT},
    {"Current 1->2 (nA)", "Total synaptic current 1->2", DefaultGUIModel::OUTPUT},

    {"Voltage 1 (V)", "Membrane potential 1", DefaultGUIModel::INPUT},
    {"Voltage 2 (V)", "Membrane potential 2", DefaultGUIModel::INPUT},
    {"Scale 2->1", "Dynamic amplitude scale 2->1", DefaultGUIModel::INPUT},
    {"Offset 2->1", "Dynamic amplitude offset 2->1", DefaultGUIModel::INPUT},
    {"Scale 1->2", "Dynamic amplitude scale 1->2", DefaultGUIModel::INPUT},
    {"Offset 1->2", "Dynamic amplitude offset 1->2", DefaultGUIModel::INPUT},
    {"Max 1 (V)", "Dynamic max 1", DefaultGUIModel::INPUT},
    {"Min 1 (V)", "Dynamic min 1", DefaultGUIModel::INPUT},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

BidirectionalChemicalSynapseGenetic::BidirectionalChemicalSynapseGenetic(void)
    : DefaultGUIModel("RTHybrid Bidirectional Chemical Synapse Genetic", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>RTHybrid Bidirectional Chemical Synapse Genetic</b></p>");
  DefaultGUIModel::createGUI(vars, num_vars);
  initParameters();
  customizeGUI();
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

double BidirectionalChemicalSynapseGenetic::sm_chemical_synapse_m(double m_slow, double v_pre, const ChemicalSynapseParams &params)
{
  return ((params.k1 * (1.0 - m_slow)) /
          (1.0 + exp(params.s_slow * (params.v_slow - v_pre)))) -
         (params.k2 * m_slow);
}

double BidirectionalChemicalSynapseGenetic::compute_i_slow(double &m_slow, double v_pre, double v_post, const ChemicalSynapseParams &params)
{
  runge_kutta_65(sm_chemical_synapse_m, m_slow, v_pre, dt, params);
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

  if (dynamic_min_max_1)
  {
    max_1 = input(8);
    min_1 = input(9);
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

  const bool aux_RT_storing = RT_storing.load(std::memory_order_acquire);
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
    if (dynamic_offsets)
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
    if (dynamic_offsets)
      new_offset = input(5) * 1000.0;

    if (new_scale != curr_scale || new_offset != curr_offset)
      scaling_factors_12_idx.store(new_scaling_factors_idx, std::memory_order_release);

    v1_scaled = v1 * new_scale + new_offset;

    const ChemicalSynapseParams &curr_params_12 = params_12[curr_synapse_idx];
    if (use_i_slow_12)
      val_i_slow_12 = compute_i_slow(m_slow_12[curr_synapse_idx], v1_scaled, v2, curr_params_12);
    if (use_i_fast_12)
      val_i_fast_12 = compute_i_fast(v1_scaled, v2, curr_params_12);
  }

  if (aux_RT_storing)
  {
    if (storing_idx < num_elements)
    {
      if (use_syn_21)
      {
        v2_scaled_sig[storing_idx] = v2_scaled;
        if (use_i_fast_21)
          i_fast_sig_21[storing_idx] = val_i_fast_21;
        if (use_i_slow_21)
          i_slow_21[storing_idx] = val_i_slow_21;
      }
      if (use_syn_12)
      {
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
  num_generations = 30u;
  population_size = 30u;

  search_phase = 1u;
  i_max_21 = 10.0;
  i_min_21 = -10.0;
  i_max_12 = 10.0;
  i_min_12 = -10.0;
  i_ranges_from_neuron = 1u;
  dt = 0.001;

  max_1 = 0.0;
  min_1 = 0.0;
  dynamic_min_max_1 = 0u;
  dynamic_offsets = 0u;
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
  genetic_running = false;

  synapse_idx.store(0, std::memory_order_relaxed);
  generations_completed = 0.0;
  individuals_completed = 0.0;
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
  {
    period = RT::System::getInstance()->getPeriod() * 1e-6; // ms

    setState("Genetic generations completed", generations_completed);
    setState("Genetic individuals of the generation completed", individuals_completed);

    setParameter("Individual evaluation time genetic (s)", evaluation_time);
    setParameter("Individual stabilization time genetic (s)", stabilization_time);
    setParameter("Genetic num generations", num_generations);
    setParameter("Genetic population size", population_size);
    setParameter("dt", dt);

    setParameter("Search phase genetic (1/0)", search_phase);
    setParameter("Current ranges to achieve are from scale of neuron (1/2)", i_ranges_from_neuron);
    setParameter("Current max to achieve genetic 2->1", i_max_21);
    setParameter("Current min to achieve genetic 2->1", i_min_21);
    setParameter("Current max to achieve genetic 1->2", i_max_12);
    setParameter("Current min to achieve genetic 1->2", i_min_12);

    setParameter("Dynamic min and max 1 (1/0)", dynamic_min_max_1);
    setParameter("Max 1 (V)", max_1);
    setParameter("Min 1 (V)", min_1);

    setParameter("Scale 2->1", scale_21_gui);
    setParameter("Offset 2->1", offset_21[scaling_factors_21_idx.load(std::memory_order_relaxed)] / 1000.0);
    setParameter("Scale 1->2", scale_12_gui);
    setParameter("Offset 1->2", offset_12[scaling_factors_12_idx.load(std::memory_order_relaxed)] / 1000.0);

    setParameter("Dynamic offsets (1/0)", dynamic_offsets);

    update_params_gui();

    setParameter("Use I_fast 2->1 (1/0)", use_i_fast_21);
    setParameter("Use I_slow 2->1 (1/0)", use_i_slow_21);
    setParameter("Use I_fast 1->2 (1/0)", use_i_fast_12);
    setParameter("Use I_slow 1->2 (1/0)", use_i_slow_12);

    break;
  }
  case MODIFY:
  {
    dt = getParameter("dt").toDouble();

    size_t curr_scaling_factors_idx, new_scaling_factors_idx;
    double curr_scale, curr_offset;

    dynamic_offsets = getParameter("Dynamic offsets (1/0)").toUInt();

    scale_21_gui = getParameter("Scale 2->1").toDouble();
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
      if (!dynamic_offsets)
        new_offset = getParameter("Offset 2->1").toDouble() * 1000.0;

      if (new_scale != curr_scale || new_offset != curr_offset)
        scaling_factors_21_idx.store(new_scaling_factors_idx, std::memory_order_release);
    }

    scale_12_gui = getParameter("Scale 1->2").toDouble();
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
      if (!dynamic_offsets)
        new_offset = getParameter("Offset 1->2").toDouble() * 1000.0;

      if (new_scale != curr_scale || new_offset != curr_offset)
        scaling_factors_12_idx.store(new_scaling_factors_idx, std::memory_order_release);
    }

    if (!genetic_running)
    {
      evaluation_time = getParameter("Individual evaluation time genetic (s)").toDouble();
      stabilization_time = getParameter("Individual stabilization time genetic (s)").toDouble();
      num_generations = getParameter("Genetic num generations").toUInt();
      population_size = getParameter("Genetic population size").toUInt();

      dynamic_min_max_1 = getParameter("Dynamic min and max 1 (1/0)").toUInt();
      if (!dynamic_min_max_1)
      {
        max_1 = getParameter("Max 1 (V)").toDouble();
        min_1 = getParameter("Min 1 (V)").toDouble();
      }

      search_phase = getParameter("Search phase genetic (1/0)").toUInt();
      i_ranges_from_neuron = getParameter("Current ranges to achieve are from scale of neuron (1/2)").toUInt();
      i_max_21 = getParameter("Current max to achieve genetic 2->1").toDouble();
      i_min_21 = getParameter("Current min to achieve genetic 2->1").toDouble();
      i_max_12 = getParameter("Current max to achieve genetic 1->2").toDouble();
      i_min_12 = getParameter("Current min to achieve genetic 1->2").toDouble();

      const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);

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
    const double new_period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
    if (new_period != period)
    {
      period = new_period;
      if (genetic_running)
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
  if (!genetic_running)
  {
    // Lógica para EMPEZAR
    if (genetic_NRT_thread.joinable())
    {
      genetic_NRT_thread.join();
      stop_genetic.store(false, std::memory_order_relaxed);
    }
    genetic_running = true;
    gentic_button->setText("Stop Genetic");
    set_params_read_only(true);
    size_t curr_scaling_factors_12_idx = scaling_factors_12_idx.load(std::memory_order_relaxed);
    genetic_NRT_thread = std::thread(&BidirectionalChemicalSynapseGenetic::NRT_genetic, this, 1.0 / period, scale_12[curr_scaling_factors_12_idx], offset_12[curr_scaling_factors_12_idx], max_1, min_1);
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
  genetic_running = false;
}

void BidirectionalChemicalSynapseGenetic::set_params_read_only(bool read_only)
{
  parameter["Individual evaluation time genetic (s)"].edit->setReadOnly(read_only);
  parameter["Individual stabilization time genetic (s)"].edit->setReadOnly(read_only);
  parameter["Genetic num generations"].edit->setReadOnly(read_only);
  parameter["Genetic population size"].edit->setReadOnly(read_only);
  parameter["Dynamic min and max 1 (1/0)"].edit->setReadOnly(read_only);
  parameter["Max 1 (V)"].edit->setReadOnly(read_only);
  parameter["Min 1 (V)"].edit->setReadOnly(read_only);

  parameter["Search phase genetic (1/0)"].edit->setReadOnly(read_only);
  parameter["Current ranges to achieve are from scale of neuron (1/2)"].edit->setReadOnly(read_only);
  parameter["Current max to achieve genetic 2->1"].edit->setReadOnly(read_only);
  parameter["Current min to achieve genetic 2->1"].edit->setReadOnly(read_only);
  parameter["Current max to achieve genetic 1->2"].edit->setReadOnly(read_only);
  parameter["Current min to achieve genetic 1->2"].edit->setReadOnly(read_only);

  parameter["E_syn 2->1"].edit->setReadOnly(read_only);
  parameter["g_fast 2->1"].edit->setReadOnly(read_only);
  parameter["s_fast 2->1"].edit->setReadOnly(read_only);
  parameter["V_fast 2->1"].edit->setReadOnly(read_only);
  parameter["g_slow 2->1"].edit->setReadOnly(read_only);
  parameter["k1 2->1"].edit->setReadOnly(read_only);
  parameter["k2 2->1"].edit->setReadOnly(read_only);
  parameter["s_slow 2->1"].edit->setReadOnly(read_only);
  parameter["V_slow 2->1"].edit->setReadOnly(read_only);
  parameter["Use I_fast 2->1 (1/0)"].edit->setReadOnly(read_only);
  parameter["Use I_slow 2->1 (1/0)"].edit->setReadOnly(read_only);

  parameter["E_syn 1->2"].edit->setReadOnly(read_only);
  parameter["g_fast 1->2"].edit->setReadOnly(read_only);
  parameter["s_fast 1->2"].edit->setReadOnly(read_only);
  parameter["V_fast 1->2"].edit->setReadOnly(read_only);
  parameter["g_slow 1->2"].edit->setReadOnly(read_only);
  parameter["k1 1->2"].edit->setReadOnly(read_only);
  parameter["k2 1->2"].edit->setReadOnly(read_only);
  parameter["s_slow 1->2"].edit->setReadOnly(read_only);
  parameter["V_slow 1->2"].edit->setReadOnly(read_only);
  parameter["Use I_fast 1->2 (1/0)"].edit->setReadOnly(read_only);
  parameter["Use I_slow 1->2 (1/0)"].edit->setReadOnly(read_only);
}

void BidirectionalChemicalSynapseGenetic::update_params_gui(void)
{
  const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_acquire);

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

void BidirectionalChemicalSynapseGenetic::set_generations_completed(double generations)
{
  generations_completed = generations;
}

void BidirectionalChemicalSynapseGenetic::set_individuals_completed(double individuals)
{
  individuals_completed = individuals;
}