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

/*
 * This is a template implementation file for a user module derived from
 * DefaultGUIModel with a custom GUI.
 */

#include "bidirectional_chemical_synapse_genetic.h"
#include <iostream>
#include <cmath>
#include <main_window.h>

enum SynapseParam
{
  SP_ESYN = 0,
  SP_G_FAST,
  SP_S_FAST,
  SP_V_FAST,
  SP_G_SLOW,
  SP_K1,
  SP_K2,
  SP_S_SLOW,
  SP_V_SLOW,
  SP_USE_I_FAST,
  SP_USE_I_SLOW,
  SP_COUNT
};

extern "C" Plugin::Object *
createRTXIPlugin(void)
{
  return new BidirectionalChemicalSynapseGenetic();
}

static DefaultGUIModel::variable_t vars[] = {
    {"Burst duration (s)", "-1 to use dynamic input", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Scale 2->1", "-1 to use dynamic input", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Offset 2->1", "-1 to use dynamic input", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Scale 1->2", "-1 to use dynamic input", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Offset 1->2", "-1 to use dynamic input", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Dynamic scaling (1/0)", "1=dynamic (input), 0=static (GUI)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

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
    {"Use I_fast 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Use I_slow 2->1 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

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
    {"Use I_fast 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},
    {"Use I_slow 1->2 (1/0)", "1 = Enable, 0 = Disable", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE},

    {"Current 2->1 (nA)", "Total synaptic current 2->1", DefaultGUIModel::OUTPUT},
    {"Current 1->2 (nA)", "Total synaptic current 1->2", DefaultGUIModel::OUTPUT},

    {"Voltage 1 (V)", "Membrane potential 1", DefaultGUIModel::INPUT},
    {"Voltage 2 (V)", "Membrane potential 2", DefaultGUIModel::INPUT},
    {"Scale 2->1", "Dynamic amplitude scale 2->1", DefaultGUIModel::INPUT},
    {"Offset 2->1", "Dynamic amplitude offset 2->1", DefaultGUIModel::INPUT},
    {"Scale 1->2", "Dynamic amplitude scale 1->2", DefaultGUIModel::INPUT},
    {"Offset 1->2", "Dynamic amplitude offset 1->2", DefaultGUIModel::INPUT},
    {"Burst duration (s)", "Dynamic burst duration", DefaultGUIModel::INPUT},

    {"m_slow 2->1", "Slow gating variable 2->1", DefaultGUIModel::STATE},
    {"m_slow 1->2", "Slow gating variable 1->2", DefaultGUIModel::STATE},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

BidirectionalChemicalSynapseGenetic::BidirectionalChemicalSynapseGenetic(void)
    : DefaultGUIModel("RTHybrid Bidirectional Chemical Synapse Genetic", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>RTHybrid Bidirectional Chemical Synapse Genetic</b></p>");
  DefaultGUIModel::createGUI(vars, num_vars);
  initParameters();
  update(INIT);
  refresh();
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

BidirectionalChemicalSynapseGenetic::~BidirectionalChemicalSynapseGenetic(void) {}

void BidirectionalChemicalSynapseGenetic::runge_kutta_65(double (*f)(double, double, double, double *), double &m_slow, double v_pre, double dt, double *params)
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

void BidirectionalChemicalSynapseGenetic::select_dt_neuron_model(double *dts, double *pts, unsigned int length, double pts_live, double *dt, double *pts_burst)
{
  double aux = pts_live;
  double factor = 1;
  double intpart, fractpart;
  int flag = 0;
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
  int length = 144;
  double pts_match = sec_per_burst * freq;
  double pts_burst, dt;

  double dts[] = {0.000500, 0.000600, 0.000700, 0.000800, 0.000900, 0.001000, 0.001100, 0.001200, 0.001300, 0.001400, 0.001500, 0.001600, 0.001800, 0.002000, 0.002200, 0.002500, 0.002800, 0.002900, 0.003000, 0.003100, 0.003200, 0.003300, 0.003400, 0.003500, 0.003600, 0.003700, 0.003800, 0.003900, 0.004000, 0.004100, 0.004200, 0.004300, 0.004400, 0.004500, 0.004600, 0.004700, 0.004800, 0.004900, 0.005000, 0.005100, 0.005200, 0.005400, 0.005600, 0.005800, 0.006000, 0.006200, 0.006400, 0.006600, 0.006800, 0.007000, 0.007200, 0.007400, 0.007700, 0.008000, 0.008300, 0.008600, 0.008900, 0.009200, 0.009600, 0.010000, 0.010400, 0.010900, 0.011400, 0.011900, 0.012500, 0.013100, 0.013800, 0.014600, 0.015400, 0.016300, 0.017300, 0.018500, 0.019900, 0.021500, 0.023300, 0.025500, 0.028100, 0.028400, 0.028700, 0.029000, 0.029400, 0.029800, 0.030200, 0.030600, 0.031000, 0.031400, 0.031800, 0.032200, 0.032600, 0.033000, 0.033400, 0.033900, 0.034400, 0.034900, 0.035400, 0.035900, 0.036400, 0.036900, 0.037400, 0.038000, 0.038600, 0.039200, 0.039800, 0.040400, 0.041000, 0.041700, 0.042400, 0.043100, 0.043800, 0.044500, 0.045300, 0.046100, 0.046900, 0.047700, 0.048600, 0.049500, 0.050400, 0.051400, 0.052400, 0.053400, 0.054500, 0.055600, 0.056800, 0.058000, 0.059300, 0.060600, 0.062000, 0.063400, 0.064900, 0.066500, 0.068200, 0.069900, 0.071700, 0.073600, 0.075600, 0.077700, 0.079900, 0.082300, 0.084800, 0.087500, 0.090300, 0.093300, 0.096500, 0.100000};
  double pts[] = {577638.000000, 481366.000000, 412599.000000, 357615.500000, 317880.000000, 286092.500000, 259143.333333, 237548.000000, 218869.500000, 203236.000000, 189687.000000, 177634.000000, 157897.000000, 142001.833333, 129024.142857, 113496.125000, 101304.555556, 97811.222222, 94527.400000, 91478.200000, 88619.400000, 85916.636364, 83389.636364, 81007.090909, 78743.583333, 76615.416667, 74599.250000, 72676.000000, 70859.076923, 69130.846154, 67476.642857, 65907.357143, 64402.666667, 62971.466667, 61602.533333, 60286.187500, 59030.250000, 57825.562500, 56664.411765, 55553.294118, 54485.000000, 52463.222222, 50586.263158, 48841.842105, 47211.050000, 45685.666667, 44255.818182, 42914.772727, 41650.739130, 40459.083333, 39335.208333, 38270.680000, 36778.346154, 35398.000000, 34117.571429, 32926.517241, 31815.833333, 30777.612903, 29493.939394, 28313.588235, 27223.638889, 25974.405405, 24834.410256, 23790.268293, 22647.767442, 21609.977778, 20513.166667, 19388.627451, 18381.132075, 17365.719298, 16361.600000, 15299.937500, 14223.202899, 13164.400000, 12147.123457, 11098.876404, 10071.693878, 9965.282828, 9861.100000, 9759.059406, 9626.242718, 9497.009615, 9371.179245, 9248.672897, 9129.293578, 9012.981818, 8899.594595, 8789.000000, 8681.149123, 8575.896552, 8473.170940, 8348.176471, 8226.809917, 8108.934426, 7994.379032, 7883.015873, 7774.710938, 7669.348837, 7566.801527, 7447.308271, 7331.525926, 7219.291971, 7110.435714, 7004.816901, 6902.298611, 6786.417808, 6674.355705, 6565.940397, 6460.987013, 6359.339744, 6247.012579, 6138.592593, 6033.866667, 5932.660714, 5822.777778, 5716.896552, 5614.796610, 5505.541436, 5400.467391, 5299.324468, 5192.348958, 5089.615385, 4982.075000, 4878.985294, 4772.014354, 4669.633803, 4564.178899, 4463.381166, 4360.214912, 4255.294872, 4149.216667, 4048.296748, 3946.654762, 3844.764479, 3743.041353, 3641.872263, 3541.583630, 3438.300000, 3336.926421, 3233.951299, 3133.666667, 3032.899696, 2932.320588, 2829.684659};
  select_dt_neuron_model(dts, pts, length, pts_match, &dt, &pts_burst);

  return pts_burst;
}

double BidirectionalChemicalSynapseGenetic::sm_chemical_synapse_m(double m_slow, double v_pre, double *params)
{
  return ((params[SP_K1] * (1.0 - m_slow)) /
          (1.0 + exp(params[SP_S_SLOW] * (params[SP_V_SLOW] - v_pre)))) -
         (params[SP_K2] * m_slow);
}

double BidirectionalChemicalSynapseGenetic::compute_synapse_current(double &m_slow, double v_pre, double v_post, double *params)
{
  double i_syn = 0.0;

  if (params[SP_USE_I_SLOW] > 0.5)
  {
    for (int i = 0; i < s_points; i++)
    {
      runge_kutta_65(sm_chemical_synapse_m, m_slow, v_pre, dt, params);
    }

    i_syn += params[SP_G_SLOW] * m_slow * (v_post - params[SP_ESYN]);
  }

  if (params[SP_USE_I_FAST] > 0.5)
  {
    i_syn += (params[SP_G_FAST] * (v_post - params[SP_ESYN])) /
             (1.0 + exp(params[SP_S_FAST] * (params[SP_V_FAST] - v_pre)));
  }

  return i_syn;
}

void BidirectionalChemicalSynapseGenetic::execute(void)
{
  if (burst_duration_gui <= 0.0)
  {
    double new_burst_duration = input(6);
    if (new_burst_duration != last_burst_duration)
    {
      last_burst_duration = new_burst_duration;
      s_points = (int)(set_pts_burst(last_burst_duration) / (last_burst_duration * freq));
      if (s_points < 1)
        s_points = 1;
    }
  }

  bool use_syn_21 = (params_21[SP_USE_I_FAST] > 0.5) || (params_21[SP_USE_I_SLOW] > 0.5);
  bool use_syn_12 = (params_12[SP_USE_I_FAST] > 0.5) || (params_12[SP_USE_I_SLOW] > 0.5);

  double v1 = input(0) * 1000.0;
  double v2 = input(1) * 1000.0;

  double v2_scaled;
  if (use_syn_21)
  {
    double scale_21, offset_21;
    if (dynamic_scaling > 0.5)
    {
      scale_21 = input(2);
      offset_21 = input(3) * 1000.0;
    }
    else
    {
      scale_21 = scale_21_gui;
      offset_21 = offset_21_gui;
    }

    if (scale_21 == 0.0)
    {
      scale_21 = 1.0;
      offset_21 = 0.0;
    }
    v2_scaled = v2 * scale_21 + offset_21;
  }

  double v1_scaled;
  if (use_syn_12)
  {
    double scale_12, offset_12;
    if (dynamic_scaling > 0.5)
    {
      scale_12 = input(4);
      offset_12 = input(5) * 1000.0;
    }
    else
    {
      scale_12 = scale_12_gui;
      offset_12 = offset_12_gui;
    }

    if (scale_12 == 0.0)
    {
      scale_12 = 1.0;
      offset_12 = 0.0;
    }
    v1_scaled = v1 * scale_12 + offset_12;
  }

  output(0) = compute_synapse_current(m_slow_21, v2_scaled, v1, params_21);
  output(1) = compute_synapse_current(m_slow_12, v1_scaled, v2, params_12);
}

void BidirectionalChemicalSynapseGenetic::initParameters(void)
{
  burst_duration_gui = 1.0;
  last_burst_duration = burst_duration_gui;

  dynamic_scaling = 0.0;
  scale_21_gui = 1.0;
  offset_21_gui = 0.0;
  scale_12_gui = 1.0;
  offset_12_gui = 0.0;

  m_slow_21 = 0.0;
  m_slow_12 = 0.0;

  params_21[SP_ESYN] = -1.92;
  params_21[SP_G_FAST] = 0.046;
  params_21[SP_S_FAST] = 0.44;
  params_21[SP_V_FAST] = -1.66;
  params_21[SP_G_SLOW] = 0.208;
  params_21[SP_K1] = 0.74;
  params_21[SP_K2] = 0.007;
  params_21[SP_S_SLOW] = 1.0;
  params_21[SP_V_SLOW] = -1.74;
  params_21[SP_USE_I_FAST] = 1.0;
  params_21[SP_USE_I_SLOW] = 1.0;

  params_12[SP_ESYN] = -1.92;
  params_12[SP_G_FAST] = 0.046;
  params_12[SP_S_FAST] = 0.44;
  params_12[SP_V_FAST] = -1.66;
  params_12[SP_G_SLOW] = 0.208;
  params_12[SP_K1] = 0.74;
  params_12[SP_K2] = 0.007;
  params_12[SP_S_SLOW] = 1.0;
  params_12[SP_V_SLOW] = -1.74;
  params_12[SP_USE_I_FAST] = 1.0;
  params_12[SP_USE_I_SLOW] = 1.0;
}

void BidirectionalChemicalSynapseGenetic::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag)
  {
  case INIT:
    period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
    freq = 1.0 / (period * 1e-3);
    s_points = (int)(set_pts_burst(last_burst_duration) / (last_burst_duration * freq));
    if (s_points == 0)
      s_points = 1;

    setParameter("Burst duration (s)", burst_duration_gui);

    setParameter("Dynamic scaling (1/0)", dynamic_scaling);
    setParameter("Scale 2->1", scale_21_gui);
    setParameter("Offset 2->1", offset_21_gui);
    setParameter("Scale 1->2", scale_12_gui);
    setParameter("Offset 1->2", offset_12_gui);

    setParameter("E_syn 2->1", params_21[SP_ESYN]);
    setParameter("g_fast 2->1", params_21[SP_G_FAST]);
    setParameter("s_fast 2->1", params_21[SP_S_FAST]);
    setParameter("V_fast 2->1", params_21[SP_V_FAST]);
    setParameter("g_slow 2->1", params_21[SP_G_SLOW]);
    setParameter("k1 2->1", params_21[SP_K1]);
    setParameter("k2 2->1", params_21[SP_K2]);
    setParameter("s_slow 2->1", params_21[SP_S_SLOW]);
    setParameter("V_slow 2->1", params_21[SP_V_SLOW]);
    setParameter("Use I_fast 2->1 (1/0)", params_21[SP_USE_I_FAST]);
    setParameter("Use I_slow 2->1 (1/0)", params_21[SP_USE_I_SLOW]);
    setState("m_slow 2->1", m_slow_21);

    setParameter("E_syn 1->2", params_12[SP_ESYN]);
    setParameter("g_fast 1->2", params_12[SP_G_FAST]);
    setParameter("s_fast 1->2", params_12[SP_S_FAST]);
    setParameter("V_fast 1->2", params_12[SP_V_FAST]);
    setParameter("g_slow 1->2", params_12[SP_G_SLOW]);
    setParameter("k1 1->2", params_12[SP_K1]);
    setParameter("k2 1->2", params_12[SP_K2]);
    setParameter("s_slow 1->2", params_12[SP_S_SLOW]);
    setParameter("V_slow 1->2", params_12[SP_V_SLOW]);
    setParameter("Use I_fast 1->2 (1/0)", params_12[SP_USE_I_FAST]);
    setParameter("Use I_slow 1->2 (1/0)", params_12[SP_USE_I_SLOW]);
    setState("m_slow 1->2", m_slow_12);

    break;

  case MODIFY:
    burst_duration_gui = getParameter("Burst duration (s)").toDouble();
    if ((burst_duration_gui > 0.0) && (burst_duration_gui != last_burst_duration))
    {
      last_burst_duration = burst_duration_gui;
      s_points = (int)(set_pts_burst(last_burst_duration) / (last_burst_duration * freq));
      if (s_points < 1)
        s_points = 1;
    }

    dynamic_scaling = getParameter("Dynamic scaling (1/0)").toDouble();
    scale_21_gui = getParameter("Scale 2->1").toDouble();
    offset_21_gui = getParameter("Offset 2->1").toDouble() * 1000.0;
    scale_12_gui = getParameter("Scale 1->2").toDouble();
    offset_12_gui = getParameter("Offset 1->2").toDouble() * 1000.0;

    params_21[SP_ESYN] = getParameter("E_syn 2->1").toDouble();
    params_21[SP_G_FAST] = getParameter("g_fast 2->1").toDouble();
    params_21[SP_S_FAST] = getParameter("s_fast 2->1").toDouble();
    params_21[SP_V_FAST] = getParameter("V_fast 2->1").toDouble();
    params_21[SP_G_SLOW] = getParameter("g_slow 2->1").toDouble();
    params_21[SP_K1] = getParameter("k1 2->1").toDouble();
    params_21[SP_K2] = getParameter("k2 2->1").toDouble();
    params_21[SP_S_SLOW] = getParameter("s_slow 2->1").toDouble();
    params_21[SP_V_SLOW] = getParameter("V_slow 2->1").toDouble();
    params_21[SP_USE_I_FAST] = getParameter("Use I_fast 2->1 (1/0)").toDouble();
    params_21[SP_USE_I_SLOW] = getParameter("Use I_slow 2->1 (1/0)").toDouble();

    params_12[SP_ESYN] = getParameter("E_syn 1->2").toDouble();
    params_12[SP_G_FAST] = getParameter("g_fast 1->2").toDouble();
    params_12[SP_S_FAST] = getParameter("s_fast 1->2").toDouble();
    params_12[SP_V_FAST] = getParameter("V_fast 1->2").toDouble();
    params_12[SP_G_SLOW] = getParameter("g_slow 1->2").toDouble();
    params_12[SP_K1] = getParameter("k1 1->2").toDouble();
    params_12[SP_K2] = getParameter("k2 1->2").toDouble();
    params_12[SP_S_SLOW] = getParameter("s_slow 1->2").toDouble();
    params_12[SP_V_SLOW] = getParameter("V_slow 1->2").toDouble();
    params_12[SP_USE_I_FAST] = getParameter("Use I_fast 1->2 (1/0)").toDouble();
    params_12[SP_USE_I_SLOW] = getParameter("Use I_slow 1->2 (1/0)").toDouble();

    break;

  case UNPAUSE:
    break;

  case PERIOD:
    period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
    freq = 1.0 / (period * 1e-3);
    s_points = (int)(set_pts_burst(last_burst_duration) / (last_burst_duration * freq));
    if (s_points == 0)
      s_points = 1;
    break;

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
  QGroupBox *button_group = new QGroupBox;

  QPushButton *abutton = new QPushButton("Button A");
  QPushButton *bbutton = new QPushButton("Button B");
  QHBoxLayout *button_layout = new QHBoxLayout;
  button_group->setLayout(button_layout);
  button_layout->addWidget(abutton);
  button_layout->addWidget(bbutton);
  QObject::connect(abutton, SIGNAL(clicked()), this, SLOT(aBttn_event()));
  QObject::connect(bbutton, SIGNAL(clicked()), this, SLOT(bBttn_event()));

  customlayout->addWidget(button_group, 0, 0);
  setLayout(customlayout);
}

void BidirectionalChemicalSynapseGenetic::aBttn_event(void) {}
void BidirectionalChemicalSynapseGenetic::bBttn_event(void) {}