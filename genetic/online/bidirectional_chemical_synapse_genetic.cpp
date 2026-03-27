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
    // Genetic basic parameters
    {"Genetic generations completed", "", DefaultGUIModel::STATE},
    {"Genetic individuals of the generation completed", "", DefaultGUIModel::STATE},
    {"Genetic num generations", "Number of generations for the genetic algorithm", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Genetic population size", "Population size for the genetic algorithm", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Genetic individual evaluation time (s)", "Individual evaluation time; does not include stabilization time", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Genetic individual stabilization time (s)", "Stabilization time for each individual; not included in evaluation time", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Genetic search phase (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Genetic current max to achieve 1->2", "In the scale of 2", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Genetic current min to achieve 1->2", "In the scale of 2", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Genetic current max to achieve 2->1", "In the scale of 1", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Genetic current min to achieve 2->1", "In the scale of 1", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    // Genetic aux parameters
    {"Dynamic voltage min and max 1 (1/0)", "1 = Enable, 0 = Disable; necessary for the genetic", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Voltage max 1", "Necessary for the genetic", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Voltage min 1", "Necessary for the genetic", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Dynamic voltage min and max 2 (1/0)", "1 = Enable, 0 = Disable; necessary for the genetic", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Voltage max 2", "Necessary for the genetic", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Voltage min 2", "Necessary for the genetic", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    // Integration timestep
    {"Step", "Integration timestep for synapse model", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    // Use fast/slow currents
    {"Use I_fast 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_slow 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
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

    {"Current 1->2", "Total synaptic current 1->2", DefaultGUIModel::OUTPUT},
    {"Current 2->1", "Total synaptic current 2->1", DefaultGUIModel::OUTPUT},

    {"Voltage 1", "Membrane potential 1", DefaultGUIModel::INPUT},
    {"Voltage 2", "Membrane potential 2", DefaultGUIModel::INPUT},
    {"Voltage max 1", "Dynamic max 1", DefaultGUIModel::INPUT},
    {"Voltage min 1", "Dynamic min 1", DefaultGUIModel::INPUT},
    {"Voltage max 2", "Dynamic max 2", DefaultGUIModel::INPUT},
    {"Voltage min 2", "Dynamic min 2", DefaultGUIModel::INPUT},
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

void BidirectionalChemicalSynapseGenetic::runge_kutta_65(double (*f)(double, double, const ChemicalSynapseParams &), double &m_slow, double v_pre, double step, const ChemicalSynapseParams &params)
{
  double apoyo, retorno;
  double k[6];

  retorno = (*f)(m_slow, v_pre, params);
  k[0] = step * retorno;
  apoyo = m_slow + k[0] * 0.2;

  retorno = (*f)(apoyo, v_pre, params);
  k[1] = step * retorno;
  apoyo = m_slow + k[0] * 0.075 + k[1] * 0.225;

  retorno = (*f)(apoyo, v_pre, params);
  k[2] = step * retorno;
  apoyo = m_slow + k[0] * 0.3 - k[1] * 0.9 + k[2] * 1.2;

  retorno = (*f)(apoyo, v_pre, params);
  k[3] = step * retorno;
  apoyo = m_slow + k[0] * 0.075 + k[1] * 0.675 - k[2] * 0.6 + k[3] * 0.75;

  retorno = (*f)(apoyo, v_pre, params);
  k[4] = step * retorno;
  apoyo = m_slow + k[0] * 0.660493827160493 + k[1] * 2.5 - k[2] * 5.185185185185185 + k[3] * 3.888888888888889 - k[4] * 0.864197530864197;

  retorno = (*f)(apoyo, v_pre, params);
  k[5] = step * retorno;

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
  runge_kutta_65(sm_chemical_synapse_m, m_slow, v_pre, step, params);
  return params.g_slow * m_slow * (v_post - params.e_syn);
}

double BidirectionalChemicalSynapseGenetic::compute_i_fast(double v_pre, double v_post, const ChemicalSynapseParams &params)
{
  return (params.g_fast * (v_post - params.e_syn)) /
         (1.0 + exp(params.s_fast * (params.v_fast - v_pre)));
}

void BidirectionalChemicalSynapseGenetic::execute(void)
{
  const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
  const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

  if (!genetic_running)
  {
    if (dynamic_v_min_max_1)
    {
      v_max_1 = input(2);
      v_min_1 = input(3);
    }

    if (dynamic_v_min_max_2)
    {
      v_max_2 = input(4);
      v_min_2 = input(5);
    }
  }

  double v1, v2;
  if (use_syn_12 || use_syn_21)
  {
    v1 = input(0);
    v2 = input(1);
  }

  double val_i_slow_12 = 0.0, val_i_fast_12 = 0.0, val_i_slow_21 = 0.0, val_i_fast_21 = 0.0;

  const bool aux_RT_storing = RT_storing.load(std::memory_order_acquire);
  const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_acquire);

  if (use_syn_12)
  {
    const ChemicalSynapseParams &curr_params_12 = params_12[curr_synapse_idx];
    if (use_i_slow_12)
      val_i_slow_12 = compute_i_slow(m_slow_12[curr_synapse_idx], v1, v2, curr_params_12);
    if (use_i_fast_12)
      val_i_fast_12 = compute_i_fast(v1, v2, curr_params_12);
  }
  output(0) = val_i_fast_12 + val_i_slow_12;

  if (use_syn_21)
  {
    const ChemicalSynapseParams &curr_params_21 = params_21[curr_synapse_idx];
    if (use_i_slow_21)
      val_i_slow_21 = compute_i_slow(m_slow_21[curr_synapse_idx], v2, v1, curr_params_21);
    if (use_i_fast_21)
      val_i_fast_21 = compute_i_fast(v2, v1, curr_params_21);
  }
  output(1) = val_i_fast_21 + val_i_slow_21;

  if (aux_RT_storing)
  {
    if (storing_idx < num_elements)
    {
      if (use_syn_12)
      {
        v_sig_1[storing_idx] = v1;
        if (use_i_fast_12)
          i_fast_sig_12[storing_idx] = val_i_fast_12;
        if (use_i_slow_12)
          i_slow_sig_12[storing_idx] = val_i_slow_12;
      }
      if (use_syn_21)
      {
        v_sig_2[storing_idx] = v2;
        if (use_i_fast_21)
          i_fast_sig_21[storing_idx] = val_i_fast_21;
        if (use_i_slow_21)
          i_slow_sig_21[storing_idx] = val_i_slow_21;
      }
      storing_idx++;
    }
    else
    {
      RT_storing.store(false, std::memory_order_release);
    }
  }
}

void BidirectionalChemicalSynapseGenetic::initParameters(void)
{
  generations_completed = 0.0;
  individuals_completed = 0.0;
  num_generations = 30u;
  population_size = 30u;
  evaluation_time = 2.0;
  stabilization_time = 0.2;
  search_phase = 1u;
  i_max_12 = 2.0;
  i_min_12 = -2.0;
  i_max_21 = 2.0;
  i_min_21 = -2.0;
  dynamic_v_min_max_1 = 0u;
  v_max_1 = 2.0;
  v_min_1 = -2.0;
  dynamic_v_min_max_2 = 0u;
  v_max_2 = 2.0;
  v_min_2 = -2.0;

  step = 0.001;

  use_i_fast_12 = 1u;
  use_i_slow_12 = 1u;
  use_i_fast_21 = 1u;
  use_i_slow_21 = 1u;

  init_syn_params_and_vars(params_12[0], m_slow_12[0]);
  init_syn_params_and_vars(params_21[0], m_slow_21[0]);

  synapse_idx.store(0, std::memory_order_relaxed);

  stop_genetic.store(false, std::memory_order_relaxed);
  RT_storing.store(false, std::memory_order_relaxed);
  genetic_running = false;
}

void BidirectionalChemicalSynapseGenetic::init_syn_params_and_vars(ChemicalSynapseParams &params, double &m_slow)
{
  params.e_syn = -1.92;
  params.g_fast = 0.046;
  params.s_fast = 0.44;
  params.v_fast = 0.0;
  params.g_slow = 0.208;
  params.k1 = 0.7;
  params.k2 = 0.7;
  params.s_slow = 1.0;
  params.v_slow = 0.0;

  m_slow = 0.0;
}

void BidirectionalChemicalSynapseGenetic::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag)
  {
  case INIT:
  {
    period = RT::System::getInstance()->getPeriod() * 1e-9; // s

    setState("Genetic generations completed", generations_completed);
    setState("Genetic individuals of the generation completed", individuals_completed);

    setParameter("Genetic num generations", num_generations);
    setParameter("Genetic population size", population_size);
    setParameter("Genetic individual evaluation time (s)", evaluation_time);
    setParameter("Genetic individual stabilization time (s)", stabilization_time);
    setParameter("Genetic search phase (1/0)", search_phase);
    setParameter("Genetic current max to achieve 1->2", i_max_12);
    setParameter("Genetic current min to achieve 1->2", i_min_12);
    setParameter("Genetic current max to achieve 2->1", i_max_21);
    setParameter("Genetic current min to achieve 2->1", i_min_21);

    setParameter("Dynamic voltage min and max 1 (1/0)", dynamic_v_min_max_1);
    setParameter("Voltage max 1", v_max_1);
    setParameter("Voltage min 1", v_min_1);
    setParameter("Dynamic voltage min and max 2 (1/0)", dynamic_v_min_max_2);
    setParameter("Voltage max 2", v_max_2);
    setParameter("Voltage min 2", v_min_2);

    setParameter("Step", step);

    setParameter("Use I_fast 1->2 (1/0)", use_i_fast_12);
    setParameter("Use I_slow 1->2 (1/0)", use_i_slow_12);
    setParameter("Use I_fast 2->1 (1/0)", use_i_fast_21);
    setParameter("Use I_slow 2->1 (1/0)", use_i_slow_21);

    update_params_gui();

    break;
  }
  case MODIFY:
  {
    if (!genetic_running)
    {
      num_generations = getParameter("Genetic num generations").toUInt();
      population_size = getParameter("Genetic population size").toUInt();
      evaluation_time = getParameter("Genetic individual evaluation time (s)").toDouble();
      stabilization_time = getParameter("Genetic individual stabilization time (s)").toDouble();
      search_phase = getParameter("Genetic search phase (1/0)").toUInt();
      i_max_12 = getParameter("Genetic current max to achieve 1->2").toDouble();
      i_min_12 = getParameter("Genetic current min to achieve 1->2").toDouble();
      i_max_21 = getParameter("Genetic current max to achieve 2->1").toDouble();
      i_min_21 = getParameter("Genetic current min to achieve 2->1").toDouble();

      dynamic_v_min_max_1 = getParameter("Dynamic voltage min and max 1 (1/0)").toUInt();
      if (!dynamic_v_min_max_1)
      {
        v_max_1 = getParameter("Voltage max 1").toDouble();
        v_min_1 = getParameter("Voltage min 1").toDouble();
      }
      dynamic_v_min_max_2 = getParameter("Dynamic voltage min and max 2 (1/0)").toUInt();
      if (!dynamic_v_min_max_2)
      {
        v_max_2 = getParameter("Voltage max 2").toDouble();
        v_min_2 = getParameter("Voltage min 2").toDouble();
      }

      step = getParameter("Step").toDouble();

      use_i_fast_12 = getParameter("Use I_fast 1->2 (1/0)").toUInt();
      use_i_slow_12 = getParameter("Use I_slow 1->2 (1/0)").toUInt();
      use_i_fast_21 = getParameter("Use I_fast 2->1 (1/0)").toUInt();
      use_i_slow_21 = getParameter("Use I_slow 2->1 (1/0)").toUInt();

      const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);

      params_12[curr_synapse_idx].e_syn = getParameter("E_syn 1->2").toDouble();
      params_12[curr_synapse_idx].g_fast = getParameter("g_fast 1->2").toDouble();
      params_12[curr_synapse_idx].s_fast = getParameter("s_fast 1->2").toDouble();
      params_12[curr_synapse_idx].v_fast = getParameter("V_fast 1->2").toDouble();
      params_12[curr_synapse_idx].g_slow = getParameter("g_slow 1->2").toDouble();
      params_12[curr_synapse_idx].k1 = getParameter("k1 1->2").toDouble();
      params_12[curr_synapse_idx].k2 = getParameter("k2 1->2").toDouble();
      params_12[curr_synapse_idx].s_slow = getParameter("s_slow 1->2").toDouble();
      params_12[curr_synapse_idx].v_slow = getParameter("V_slow 1->2").toDouble();

      params_21[curr_synapse_idx].e_syn = getParameter("E_syn 2->1").toDouble();
      params_21[curr_synapse_idx].g_fast = getParameter("g_fast 2->1").toDouble();
      params_21[curr_synapse_idx].s_fast = getParameter("s_fast 2->1").toDouble();
      params_21[curr_synapse_idx].v_fast = getParameter("V_fast 2->1").toDouble();
      params_21[curr_synapse_idx].g_slow = getParameter("g_slow 2->1").toDouble();
      params_21[curr_synapse_idx].k1 = getParameter("k1 2->1").toDouble();
      params_21[curr_synapse_idx].k2 = getParameter("k2 2->1").toDouble();
      params_21[curr_synapse_idx].s_slow = getParameter("s_slow 2->1").toDouble();
      params_21[curr_synapse_idx].v_slow = getParameter("V_slow 2->1").toDouble();
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
    genetic_NRT_thread = std::thread(&BidirectionalChemicalSynapseGenetic::NRT_genetic, this, period, v_max_1, v_min_1, v_max_2, v_min_2);
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
  for (std::pair<const QString, param_t> &elemento : parameter)
  {
    elemento.second.edit->setReadOnly(read_only);
  }
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