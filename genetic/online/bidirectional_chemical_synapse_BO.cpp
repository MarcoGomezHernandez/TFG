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

#include "bidirectional_chemical_synapse_BO.h"
#include <chrono>
#include <algorithm>
#include <main_window.h>

namespace ModulePrivateConfig
{
  // Default cutoff frequency (kHz) used when separating I_fast/I_slow for BO scoring.
  static constexpr double FILTER_FC = 0.3;

  // Small margin used when clamping the slow gating variable (avoids sticking at bounds).
  static constexpr double M_SLOW_MARGIN = 1e-6;
}

namespace ModuleConstants
{
  // Slightly expanded clamp to keep numerical integrator stable.
  static constexpr double M_SLOW_MIN = -ModulePrivateConfig::M_SLOW_MARGIN;
  static constexpr double M_SLOW_MAX = 1.0 + ModulePrivateConfig::M_SLOW_MARGIN;
}

extern "C" Plugin::Object *
createRTXIPlugin(void)
{
  // RTXI plugin factory.
  return new BidirectionalChemicalSynapseBO();
}

// GUI variables exposed to RTXI.
// This list includes:
// - BO configuration (initial samples, iterations, times, search phase, cutoffs)
// - Target current bounds to achieve (expected min/max for both directions)
// - Voltage bounds (fixed or dynamic) used for BO parameter-range initialization
// - Output clamps (min/max current) for safety
// - Component toggles (use fast/slow per direction)
// - The actual synapse parameters (E_syn, g_*, s_*, V_*, k1/k2) per direction
// - IO signals (voltages in, currents out)
static DefaultGUIModel::variable_t vars[] = {
    {"BO evaluations completed", "Finishes when this is initial samples + iterations", DefaultGUIModel::STATE},
    {"BO initial samples", "Number of initialization samples for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"BO iterations", "Number of BO iterations after initial sampling", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"BO evaluation time (ms)", "Time to record signals per evaluation", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO stabilization time (ms)", "Wait time after setting params before recording", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO search phase (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"BO current min to achieve 1->2 (nA)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO current max to achieve 1->2 (nA)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO current min to achieve 2->1 (nA)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO current max to achieve 2->1 (nA)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO cutoff frequency 1 (kHz)", "To separate the I_fast and I_slow for BO in synapse 1->2", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"BO cutoff frequency 2 (kHz)", "To separate the I_fast and I_slow for BO in synapse 2->1", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Dynamic voltage min and max 1 (1/0)", "1 = Enable, 0 = Disable; necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Voltage min 1 (V)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Voltage max 1 (V)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Dynamic voltage min and max 2 (1/0)", "1 = Enable, 0 = Disable; necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Voltage min 2 (V)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Voltage max 2 (V)", "Necessary for BO", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Current min 1->2 (nA)", "Fixed output clamp min for current 1->2", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Current max 1->2 (nA)", "Fixed output clamp max for current 1->2", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Current min 2->1 (nA)", "Fixed output clamp min for current 2->1", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Current max 2->1 (nA)", "Fixed output clamp max for current 2->1", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Verbose (1/0)", "Enable/disable BO candidate evaluation logging", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},

    {"factor in dt (ms) = period (ms) * factor", "Factor for calculating dt form the period; dt in ms", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Use I_fast 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_slow 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_fast 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
    {"Use I_slow 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},

    {"E_syn 1->2 (V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_fast 1->2 (nS)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_fast 1->2 (1/V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_fast 1->2 (V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_slow 1->2 (nS)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k1 1->2 (1/ms)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k2 1->2 (1/ms)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_slow 1->2 (1/V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_slow 1->2 (V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"E_syn 2->1 (V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_fast 2->1 (nS)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_fast 2->1 (1/V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_fast 2->1 (V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"g_slow 2->1 (nS)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k1 2->1 (1/ms)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"k2 2->1 (1/ms)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"s_slow 2->1 (1/V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"V_slow 2->1 (V)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Current 1->2 (nA)", "Total synaptic current 1->2", DefaultGUIModel::OUTPUT},
    {"Current 2->1 (nA)", "Total synaptic current 2->1", DefaultGUIModel::OUTPUT},

    {"Voltage 1 (V)", "Membrane potential 1", DefaultGUIModel::INPUT},
    {"Voltage 2 (V)", "Membrane potential 2", DefaultGUIModel::INPUT},
    {"Voltage min 1 (V)", "Dynamic min 1", DefaultGUIModel::INPUT},
    {"Voltage max 1 (V)", "Dynamic max 1", DefaultGUIModel::INPUT},
    {"Voltage min 2 (V)", "Dynamic min 2", DefaultGUIModel::INPUT},
    {"Voltage max 2 (V)", "Dynamic max 2", DefaultGUIModel::INPUT},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

BidirectionalChemicalSynapseBO::BidirectionalChemicalSynapseBO(void)
    : DefaultGUIModel("RTHybrid Bidirectional Chemical Synapse BO", ::vars, ::num_vars)
{
  // Construct RTXI GUI and initialize parameters.
  setWhatsThis("<p><b>RTHybrid Bidirectional Chemical Synapse BO</b></p>");
  DefaultGUIModel::createGUI(vars, num_vars);
  initParameters();
  customizeGUI();
  update(INIT);
  refresh();
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

BidirectionalChemicalSynapseBO::~BidirectionalChemicalSynapseBO(void)
{
  if (BO_NRT_thread.joinable())
  {
    // Ensure BO thread exits before destruction.
    stop_BO.store(true, std::memory_order_relaxed);
    BO_NRT_thread.join();
  }
}

void BidirectionalChemicalSynapseBO::runge_kutta_65(double (*f)(double, double, const ChemicalSynapseParams &), double &m_slow, double v_pre, double dt, const ChemicalSynapseParams &params)
{
  // Fixed-step 6-stage Runge–Kutta integrator (coefficients are hard-coded).
  // Used to integrate the slow synaptic gating variable ODE.
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

double BidirectionalChemicalSynapseBO::sm_chemical_synapse_m(double m_slow, double v_pre, const ChemicalSynapseParams &params)
{
  // Slow gate dynamics: dm/dt = k1*(1-m)*sigmoid(...) - k2*m.
  return (params.k1 * (1.0 - m_slow) * chemical_sigmoid(params.s_slow, params.v_slow, v_pre)) -
         (params.k2 * m_slow);
}

double BidirectionalChemicalSynapseBO::compute_i_slow(double &m_slow, double v_pre, double v_post, const ChemicalSynapseParams &params)
{
  // Clamp state, integrate ODE one step, and compute slow current.
  m_slow = std::clamp(m_slow, ModuleConstants::M_SLOW_MIN, ModuleConstants::M_SLOW_MAX);
  runge_kutta_65(sm_chemical_synapse_m, m_slow, v_pre, dt, params);
  return params.g_slow * m_slow * (v_post - params.e_syn);
}

double BidirectionalChemicalSynapseBO::compute_i_fast(double v_pre, double v_post, const ChemicalSynapseParams &params)
{
  // Fast current is instantaneous (sigmoid activation).
  return params.g_fast * (v_post - params.e_syn) * chemical_sigmoid(params.s_fast, params.v_fast, v_pre);
}

void BidirectionalChemicalSynapseBO::execute(void)
{
  // Real-time loop:
  // - Reads inputs (voltages and optional dynamic bounds)
  // - Computes bidirectional synaptic currents (fast + slow)
  // - Applies output clamps
  // - Optionally records signals into pre-allocated buffers when BO requests it
  if (!BO_running)
  {
    if (dynamic_v_min_max_1)
    {
      v_min_1 = input(2);
      v_max_1 = input(3);
    }

    if (dynamic_v_min_max_2)
    {
      v_min_2 = input(4);
      v_max_2 = input(5);
    }
  }

  const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
  const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;

  double v1, v2;
  if (use_syn_12 || use_syn_21)
  {
    v1 = input(0);
    v2 = input(1);
  }

  double val_i_slow_12 = 0.0, val_i_fast_12 = 0.0, val_i_slow_21 = 0.0, val_i_fast_21 = 0.0;

  // Handshake flags from NRT -> RT.
  // Acquire ensures we see the latest `num_elements`/buffers before storing.
  const bool aux_RT_storing = RT_storing.load(std::memory_order_acquire);
  // Acquire pairs with NRT release-store so RT sees a fully written params slot.
  const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_acquire);
  // Inform NRT which slot RT is currently using (used to coordinate safe swaps).
  last_synapse_idx_read_RT.store(curr_synapse_idx, std::memory_order_relaxed);

  if (use_syn_12)
  {
    const ChemicalSynapseParams &curr_params_12 = params_12[curr_synapse_idx];
    if (use_i_slow_12)
      val_i_slow_12 = compute_i_slow(m_slow_12, v1, v2, curr_params_12);
    if (use_i_fast_12)
      val_i_fast_12 = compute_i_fast(v1, v2, curr_params_12);
  }
  const double val_i_12 = val_i_fast_12 + val_i_slow_12;
  output(0) = std::clamp(val_i_12, i_min_12, i_max_12);

  if (use_syn_21)
  {
    const ChemicalSynapseParams &curr_params_21 = params_21[curr_synapse_idx];
    if (use_i_slow_21)
      val_i_slow_21 = compute_i_slow(m_slow_21, v2, v1, curr_params_21);
    if (use_i_fast_21)
      val_i_fast_21 = compute_i_fast(v2, v1, curr_params_21);
  }
  const double val_i_21 = val_i_fast_21 + val_i_slow_21;
  output(1) = std::clamp(val_i_21, i_min_21, i_max_21);

  if (aux_RT_storing)
  {
    // Store samples until `num_elements` is reached, then clear the flag.
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
      // Release so NRT sees all buffer writes before RT_storing becomes false.
      RT_storing.store(false, std::memory_order_release);
    }
  }
}

void BidirectionalChemicalSynapseBO::initParameters(void)
{
  // Default BO configuration (modifiable via GUI when BO is stopped).
  // - initial_samples/iterations: Limbo sampling schedule
  // - evaluation_time/stabilization_time: capture timing in ms
  // - search_phase: affects scoring polarity and E_syn search region
  // - expected_i_*: target min/max currents to achieve (nA)
  // - fc_1/fc_2: cutoff (kHz) used to separate I_fast/I_slow for scoring
  evaluations_completed = 0.0;
  initial_samples = 40u;
  iterations = 200u;
  evaluation_time = 2000.0;
  stabilization_time = 1000.0;
  search_phase = 1u;
  expected_i_min_12 = 0.0;
  expected_i_max_12 = 0.0;
  expected_i_min_21 = 0.0;
  expected_i_max_21 = 0.0;
  constexpr double FILTER_FC = ModulePrivateConfig::FILTER_FC;
  fc_1 = FILTER_FC;
  fc_2 = FILTER_FC;

  // Voltage bounds used for BO range initialization (fixed or provided dynamically).
  dynamic_v_min_max_1 = 0u;
  v_min_1 = 0.0;
  v_max_1 = 0.0;

  dynamic_v_min_max_2 = 0u;
  v_min_2 = 0.0;
  v_max_2 = 0.0;

  // Output clamps (min/max) for each direction (nA).
  i_min_12 = 0.0;
  i_max_12 = 0.0;
  i_min_21 = 0.0;
  i_max_21 = 0.0;

  // Verbose flag for BO candidate evaluation logging.
  verbose.store(0u, std::memory_order_relaxed);

  // dt = period * dt_factor (ms). Changing dt affects ODE integration.
  dt_factor = 1.0;

  // Enable/disable fast and slow components per direction.
  use_i_fast_12 = 1u;
  use_i_slow_12 = 1u;
  use_i_fast_21 = 1u;
  use_i_slow_21 = 1u;

  init_syn_params_and_vars(params_12[0]);
  init_syn_params_and_vars(params_12[1]);
  init_syn_params_and_vars(params_21[0]);
  init_syn_params_and_vars(params_21[1]);

  m_slow_12 = 0.0;
  m_slow_21 = 0.0;

  // Start on slot 0; NRT swaps between 0 and 1.
  synapse_idx.store(0, std::memory_order_relaxed);
  last_synapse_idx_read_RT.store(0, std::memory_order_relaxed);

  // BO thread control + capture control.
  stop_BO.store(false, std::memory_order_relaxed);
  RT_storing.store(false, std::memory_order_relaxed);
  BO_running = false;

  // Limbo stop criterion reads this flag.
  StopFunctor::stop_BO_ptr = &stop_BO;
}

void BidirectionalChemicalSynapseBO::init_syn_params_and_vars(ChemicalSynapseParams &params)
{
  params.e_syn = 0.0;
  params.g_fast = 0.0;
  params.s_fast = 0.0;
  params.v_fast = 0.0;
  params.g_slow = 0.0;
  params.k1 = 0.0;
  params.k2 = 0.0;
  params.s_slow = 0.0;
  params.v_slow = 0.0;
}

bool BidirectionalChemicalSynapseBO::wait_until_RT_read_idx_or_stop(size_t idx_to_achieve)
{
  const std::chrono::duration<double, std::milli> active_wait_duration(BOPublicConfig::ACTIVE_WAIT_MS);

  // Active-wait until RT confirms it has read `idx_to_achieve`.
  // This prevents NRT from overwriting the slot currently being used by RT.
  while (last_synapse_idx_read_RT.load(std::memory_order_relaxed) != idx_to_achieve)
  {
    if (stop_BO.load(std::memory_order_relaxed))
      return false;

    std::this_thread::sleep_for(active_wait_duration);
  }

  return true;
}

void BidirectionalChemicalSynapseBO::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag)
  {
  case INIT:
  {
    // Initial GUI population + derived timing (period/dt).
    period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
    dt = period * dt_factor;                                // ms

    setState("BO evaluations completed", evaluations_completed);

    setParameter("BO initial samples", initial_samples);
    setParameter("BO iterations", iterations);
    setParameter("BO evaluation time (ms)", evaluation_time);
    setParameter("BO stabilization time (ms)", stabilization_time);
    setParameter("BO search phase (1/0)", search_phase);
    setParameter("BO current min to achieve 1->2 (nA)", expected_i_min_12);
    setParameter("BO current max to achieve 1->2 (nA)", expected_i_max_12);
    setParameter("BO current min to achieve 2->1 (nA)", expected_i_min_21);
    setParameter("BO current max to achieve 2->1 (nA)", expected_i_max_21);
    setParameter("BO cutoff frequency 1 (kHz)", fc_1);
    setParameter("BO cutoff frequency 2 (kHz)", fc_2);

    setParameter("Dynamic voltage min and max 1 (1/0)", dynamic_v_min_max_1);
    setParameter("Voltage min 1 (V)", v_min_1);
    setParameter("Voltage max 1 (V)", v_max_1);
    setParameter("Dynamic voltage min and max 2 (1/0)", dynamic_v_min_max_2);
    setParameter("Voltage min 2 (V)", v_min_2);
    setParameter("Voltage max 2 (V)", v_max_2);

    setParameter("Current min 1->2 (nA)", i_min_12);
    setParameter("Current max 1->2 (nA)", i_max_12);
    setParameter("Current min 2->1 (nA)", i_min_21);
    setParameter("Current max 2->1 (nA)", i_max_21);

    setParameter("Verbose (1/0)", verbose.load(std::memory_order_relaxed));

    setParameter("factor in dt (ms) = period (ms) * factor", dt_factor);

    setParameter("Use I_fast 1->2 (1/0)", use_i_fast_12);
    setParameter("Use I_slow 1->2 (1/0)", use_i_slow_12);
    setParameter("Use I_fast 2->1 (1/0)", use_i_fast_21);
    setParameter("Use I_slow 2->1 (1/0)", use_i_slow_21);

    update_params_gui();

    break;
  }
  case MODIFY:
  {
    // Apply GUI changes. Some changes are only allowed when BO is stopped to avoid conflicts with the BO thread and nonsenses i the BO.
    i_min_12 = getParameter("Current min 1->2 (nA)").toDouble();
    i_max_12 = getParameter("Current max 1->2 (nA)").toDouble();
    i_min_21 = getParameter("Current min 2->1 (nA)").toDouble();
    i_max_21 = getParameter("Current max 2->1 (nA)").toDouble();

    verbose.store(getParameter("Verbose (1/0)").toUInt(), std::memory_order_relaxed);

    if (!BO_running)
    {
      initial_samples = getParameter("BO initial samples").toUInt();
      iterations = getParameter("BO iterations").toUInt();
      evaluation_time = getParameter("BO evaluation time (ms)").toDouble();
      stabilization_time = getParameter("BO stabilization time (ms)").toDouble();
      search_phase = getParameter("BO search phase (1/0)").toUInt();
      expected_i_min_12 = getParameter("BO current min to achieve 1->2 (nA)").toDouble();
      expected_i_max_12 = getParameter("BO current max to achieve 1->2 (nA)").toDouble();
      expected_i_min_21 = getParameter("BO current min to achieve 2->1 (nA)").toDouble();
      expected_i_max_21 = getParameter("BO current max to achieve 2->1 (nA)").toDouble();
      fc_1 = getParameter("BO cutoff frequency 1 (kHz)").toDouble();
      fc_2 = getParameter("BO cutoff frequency 2 (kHz)").toDouble();

      dynamic_v_min_max_1 = getParameter("Dynamic voltage min and max 1 (1/0)").toUInt();
      const double new_v_min_1 = getParameter("Voltage min 1 (V)").toDouble();
      const double new_v_max_1 = getParameter("Voltage max 1 (V)").toDouble();
      if (!dynamic_v_min_max_1)
      {
        v_min_1 = new_v_min_1;
        v_max_1 = new_v_max_1;
      }
      dynamic_v_min_max_2 = getParameter("Dynamic voltage min and max 2 (1/0)").toUInt();
      const double new_v_min_2 = getParameter("Voltage min 2 (V)").toDouble();
      const double new_v_max_2 = getParameter("Voltage max 2 (V)").toDouble();
      if (!dynamic_v_min_max_2)
      {
        v_min_2 = new_v_min_2;
        v_max_2 = new_v_max_2;
      }

      double new_dt_factor = getParameter("factor in dt (ms) = period (ms) * factor").toDouble();
      if (new_dt_factor != dt_factor)
      {
        dt_factor = new_dt_factor;
        dt = period * dt_factor;
      }

      use_i_fast_12 = getParameter("Use I_fast 1->2 (1/0)").toUInt();
      use_i_slow_12 = getParameter("Use I_slow 1->2 (1/0)").toUInt();
      use_i_fast_21 = getParameter("Use I_fast 2->1 (1/0)").toUInt();
      use_i_slow_21 = getParameter("Use I_slow 2->1 (1/0)").toUInt();

      const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);

      params_12[curr_synapse_idx].e_syn = getParameter("E_syn 1->2 (V)").toDouble();
      params_12[curr_synapse_idx].g_fast = getParameter("g_fast 1->2 (nS)").toDouble();
      params_12[curr_synapse_idx].s_fast = getParameter("s_fast 1->2 (1/V)").toDouble();
      params_12[curr_synapse_idx].v_fast = getParameter("V_fast 1->2 (V)").toDouble();
      params_12[curr_synapse_idx].g_slow = getParameter("g_slow 1->2 (nS)").toDouble();
      params_12[curr_synapse_idx].k1 = getParameter("k1 1->2 (1/ms)").toDouble();
      params_12[curr_synapse_idx].k2 = getParameter("k2 1->2 (1/ms)").toDouble();
      params_12[curr_synapse_idx].s_slow = getParameter("s_slow 1->2 (1/V)").toDouble();
      params_12[curr_synapse_idx].v_slow = getParameter("V_slow 1->2 (V)").toDouble();

      params_21[curr_synapse_idx].e_syn = getParameter("E_syn 2->1 (V)").toDouble();
      params_21[curr_synapse_idx].g_fast = getParameter("g_fast 2->1 (nS)").toDouble();
      params_21[curr_synapse_idx].s_fast = getParameter("s_fast 2->1 (1/V)").toDouble();
      params_21[curr_synapse_idx].v_fast = getParameter("V_fast 2->1 (V)").toDouble();
      params_21[curr_synapse_idx].g_slow = getParameter("g_slow 2->1 (nS)").toDouble();
      params_21[curr_synapse_idx].k1 = getParameter("k1 2->1 (1/ms)").toDouble();
      params_21[curr_synapse_idx].k2 = getParameter("k2 2->1 (1/ms)").toDouble();
      params_21[curr_synapse_idx].s_slow = getParameter("s_slow 2->1 (1/V)").toDouble();
      params_21[curr_synapse_idx].v_slow = getParameter("V_slow 2->1 (V)").toDouble();
    }

    break;
  }

  case UNPAUSE:
    break;

  case PERIOD:
  {
    // RT period changed: dt and number of stored samples change -> stop BO.
    const double new_period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
    if (new_period != period)
    {
      period = new_period;
      dt = period * dt_factor;
      if (BO_running)
      {
        // Stop BO to keep capture length consistent.
        stop_BO.store(true, std::memory_order_relaxed); // Porque cambiaría el número de puntos a almacenar
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

void BidirectionalChemicalSynapseBO::customizeGUI(void)
{
  QGridLayout *customlayout = DefaultGUIModel::getLayout();
  // Add a simple Start/Stop button to control the BO thread.
  BO_button = new QPushButton("Start BO");
  QObject::connect(BO_button, SIGNAL(clicked()), this, SLOT(toggle_BO_event()));
  customlayout->addWidget(BO_button, 0, 0);
  setLayout(customlayout);
}

void BidirectionalChemicalSynapseBO::toggle_BO_event(void)
{
  if (!BO_running)
  {
    // Start BO: ensure any previous thread is joined and clear stop flag.
    if (BO_NRT_thread.joinable())
    {
      BO_NRT_thread.join();
      stop_BO.store(false, std::memory_order_relaxed);
    }
    BO_running = true;
    BO_button->setText("Stop BO");
    set_params_read_only(true);
    // Launch NRT optimization loop; it uses atomics to coordinate with RT.
    BO_NRT_thread = std::thread(&BidirectionalChemicalSynapseBO::NRT_BO, this, period);
  }
  else
  {
    // Request stop (actual thread exit is handled inside the BO loop).
    stop_BO.store(true, std::memory_order_relaxed);
  }
}

void BidirectionalChemicalSynapseBO::stop_BO_event_async(void)
{
  update_params_gui();
  set_params_read_only(false);
  BO_button->setText("Start BO");
  BO_running = false;
}

void BidirectionalChemicalSynapseBO::set_params_read_only(bool read_only)
{
  // Lock GUI fields while BO is running to avoid inconsistent configurations.
  parameter["BO initial samples"].edit->setReadOnly(read_only);
  parameter["BO iterations"].edit->setReadOnly(read_only);
  parameter["BO evaluation time (ms)"].edit->setReadOnly(read_only);
  parameter["BO stabilization time (ms)"].edit->setReadOnly(read_only);
  parameter["BO search phase (1/0)"].edit->setReadOnly(read_only);
  parameter["BO current min to achieve 1->2 (nA)"].edit->setReadOnly(read_only);
  parameter["BO current max to achieve 1->2 (nA)"].edit->setReadOnly(read_only);
  parameter["BO current min to achieve 2->1 (nA)"].edit->setReadOnly(read_only);
  parameter["BO current max to achieve 2->1 (nA)"].edit->setReadOnly(read_only);
  parameter["BO cutoff frequency 1 (kHz)"].edit->setReadOnly(read_only);
  parameter["BO cutoff frequency 2 (kHz)"].edit->setReadOnly(read_only);

  parameter["Dynamic voltage min and max 1 (1/0)"].edit->setReadOnly(read_only);
  parameter["Voltage min 1 (V)"].edit->setReadOnly(read_only);
  parameter["Voltage max 1 (V)"].edit->setReadOnly(read_only);
  parameter["Dynamic voltage min and max 2 (1/0)"].edit->setReadOnly(read_only);
  parameter["Voltage min 2 (V)"].edit->setReadOnly(read_only);
  parameter["Voltage max 2 (V)"].edit->setReadOnly(read_only);

  parameter["factor in dt (ms) = period (ms) * factor"].edit->setReadOnly(read_only);

  parameter["Use I_fast 1->2 (1/0)"].edit->setReadOnly(read_only);
  parameter["Use I_slow 1->2 (1/0)"].edit->setReadOnly(read_only);
  parameter["Use I_fast 2->1 (1/0)"].edit->setReadOnly(read_only);
  parameter["Use I_slow 2->1 (1/0)"].edit->setReadOnly(read_only);

  if (use_i_fast_12 || use_i_slow_12)
  {
    parameter["E_syn 1->2 (V)"].edit->setReadOnly(read_only);
    if (use_i_fast_12)
    {
      parameter["g_fast 1->2 (nS)"].edit->setReadOnly(read_only);
      parameter["s_fast 1->2 (1/V)"].edit->setReadOnly(read_only);
      parameter["V_fast 1->2 (V)"].edit->setReadOnly(read_only);
    }
    if (use_i_slow_12)
    {
      parameter["g_slow 1->2 (nS)"].edit->setReadOnly(read_only);
      parameter["k1 1->2 (1/ms)"].edit->setReadOnly(read_only);
      parameter["k2 1->2 (1/ms)"].edit->setReadOnly(read_only);
      parameter["s_slow 1->2 (1/V)"].edit->setReadOnly(read_only);
      parameter["V_slow 1->2 (V)"].edit->setReadOnly(read_only);
    }
  }

  if (use_i_fast_21 || use_i_slow_21)
  {
    parameter["E_syn 2->1 (V)"].edit->setReadOnly(read_only);
    if (use_i_fast_21)
    {
      parameter["g_fast 2->1 (nS)"].edit->setReadOnly(read_only);
      parameter["s_fast 2->1 (1/V)"].edit->setReadOnly(read_only);
      parameter["V_fast 2->1 (V)"].edit->setReadOnly(read_only);
    }
    if (use_i_slow_21)
    {
      parameter["g_slow 2->1 (nS)"].edit->setReadOnly(read_only);
      parameter["k1 2->1 (1/ms)"].edit->setReadOnly(read_only);
      parameter["k2 2->1 (1/ms)"].edit->setReadOnly(read_only);
      parameter["s_slow 2->1 (1/V)"].edit->setReadOnly(read_only);
      parameter["V_slow 2->1 (V)"].edit->setReadOnly(read_only);
    }
  }
}

void BidirectionalChemicalSynapseBO::update_params_gui(void)
{
  // Pull the currently active params slot into the GUI.
  // Acquire ensures the slot contents are fully visible after NRT publishes it.
  const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_acquire);
  // Handshake: let NRT know RT/GUI observed this slot.
  last_synapse_idx_read_RT.store(curr_synapse_idx, std::memory_order_relaxed);

  if (use_i_fast_12 || use_i_slow_12)
  {
    setParameter("E_syn 1->2 (V)", params_12[curr_synapse_idx].e_syn);
    if (use_i_fast_12)
    {
      setParameter("g_fast 1->2 (nS)", params_12[curr_synapse_idx].g_fast);
      setParameter("s_fast 1->2 (1/V)", params_12[curr_synapse_idx].s_fast);
      setParameter("V_fast 1->2 (V)", params_12[curr_synapse_idx].v_fast);
    }
    if (use_i_slow_12)
    {
      setParameter("g_slow 1->2 (nS)", params_12[curr_synapse_idx].g_slow);
      setParameter("k1 1->2 (1/ms)", params_12[curr_synapse_idx].k1);
      setParameter("k2 1->2 (1/ms)", params_12[curr_synapse_idx].k2);
      setParameter("s_slow 1->2 (1/V)", params_12[curr_synapse_idx].s_slow);
      setParameter("V_slow 1->2 (V)", params_12[curr_synapse_idx].v_slow);
    }
  }

  if (use_i_fast_21 || use_i_slow_21)
  {
    setParameter("E_syn 2->1 (V)", params_21[curr_synapse_idx].e_syn);
    if (use_i_fast_21)
    {
      setParameter("g_fast 2->1 (nS)", params_21[curr_synapse_idx].g_fast);
      setParameter("s_fast 2->1 (1/V)", params_21[curr_synapse_idx].s_fast);
      setParameter("V_fast 2->1 (V)", params_21[curr_synapse_idx].v_fast);
    }
    if (use_i_slow_21)
    {
      setParameter("g_slow 2->1 (nS)", params_21[curr_synapse_idx].g_slow);
      setParameter("k1 2->1 (1/ms)", params_21[curr_synapse_idx].k1);
      setParameter("k2 2->1 (1/ms)", params_21[curr_synapse_idx].k2);
      setParameter("s_slow 2->1 (1/V)", params_21[curr_synapse_idx].s_slow);
      setParameter("V_slow 2->1 (V)", params_21[curr_synapse_idx].v_slow);
    }
  }
}

void BidirectionalChemicalSynapseBO::set_evaluations_completed(double evals)
{
  evaluations_completed = evals;
}