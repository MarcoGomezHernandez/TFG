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

#include "rthybrid_komendantov_kononenko_1996_neuron.h"
#include <iostream>
#include <main_window.h>

#define NM_KOMENDANTOV_KONONENKO_1996_N_VARS 8
// Variables
#define NM_KOMENDANTOV_KONONENKO_1996_V 0
#define NM_KOMENDANTOV_KONONENKO_1996_M_B 1
#define NM_KOMENDANTOV_KONONENKO_1996_H_B 2
#define NM_KOMENDANTOV_KONONENKO_1996_M 3
#define NM_KOMENDANTOV_KONONENKO_1996_H 4
#define NM_KOMENDANTOV_KONONENKO_1996_N 5
#define NM_KOMENDANTOV_KONONENKO_1996_M_CA 6
#define NM_KOMENDANTOV_KONONENKO_1996_CA 7


#define NM_KOMENDANTOV_KONONENKO_1996_N_PARAMS 20
// Parameters
#define NM_KOMENDANTOV_KONONENKO_1996_DT 0
#define NM_KOMENDANTOV_KONONENKO_1996_I 1
#define NM_KOMENDANTOV_KONONENKO_1996_SYN 2
#define NM_KOMENDANTOV_KONONENKO_1996_CM 3
#define NM_KOMENDANTOV_KONONENKO_1996_G_NA_V 4
#define NM_KOMENDANTOV_KONONENKO_1996_V_NA 5
#define NM_KOMENDANTOV_KONONENKO_1996_G_K 6
#define NM_KOMENDANTOV_KONONENKO_1996_V_K 7
#define NM_KOMENDANTOV_KONONENKO_1996_G_NA 8
#define NM_KOMENDANTOV_KONONENKO_1996_G_B 9
#define NM_KOMENDANTOV_KONONENKO_1996_V_B 10
#define NM_KOMENDANTOV_KONONENKO_1996_G_NA_TTX 11
#define NM_KOMENDANTOV_KONONENKO_1996_G_K_TEA 12
#define NM_KOMENDANTOV_KONONENKO_1996_G_CA 13
#define NM_KOMENDANTOV_KONONENKO_1996_V_CA 14
#define NM_KOMENDANTOV_KONONENKO_1996_G_CA_CA 15
#define NM_KOMENDANTOV_KONONENKO_1996_K_BETA 16
#define NM_KOMENDANTOV_KONONENKO_1996_BETA 17
#define NM_KOMENDANTOV_KONONENKO_1996_RHO 18
#define NM_KOMENDANTOV_KONONENKO_1996_K_S 19

extern "C" Plugin::Object*
createRTXIPlugin(void)
{
  return new RTHybridKomendantovKononenko1996Neuron();
}

static DefaultGUIModel::variable_t vars[] = {
  {"Burst duration (s)", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"i", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"cm", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

  {"g_na_v", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"v_na", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"g_na", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"g_na_ttx", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

  {"g_k", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"g_k_tea", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"v_k", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

  {"g_b", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"v_b", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

  {"g_ca", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"v_ca", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"g_ca_ca", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"k_beta", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"beta", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"rho", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
  {"k_s", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

  { "v0 (mV)", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

  { "Vm (v)", "Membrane potential (in V)", DefaultGUIModel::OUTPUT, },
  { "Vm (mV)", "Membrane potential (in mV)", DefaultGUIModel::OUTPUT, },

  { "Isyn (nA)", "Synaptic input current (in nA)", DefaultGUIModel::INPUT,},
  { "Burst duration (s)", "", DefaultGUIModel::INPUT,},

  { "v", "", DefaultGUIModel::STATE,},
  { "s_points", "", DefaultGUIModel::STATE,},
  { "dt", "", DefaultGUIModel::STATE,},
  { "syn", "", DefaultGUIModel::STATE,},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

RTHybridKomendantovKononenko1996Neuron::RTHybridKomendantovKononenko1996Neuron(void)
  : DefaultGUIModel("RTHybrid Komendantov-Kononenko (1996) neuron model", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>RTHybrid Komendantov-Kononenko (1996) neuron model</b></p>");
  DefaultGUIModel::createGUI(vars,
                             num_vars); // this is required to create the GUI
  //customizeGUI();
  initParameters();
  update(INIT); // this is optional, you may place initialization code directly
                // into the constructor
  refresh();    // this is required to update the GUI with parameter and state
                // values
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

RTHybridKomendantovKononenko1996Neuron::~RTHybridKomendantovKononenko1996Neuron(void)
{
}

void RTHybridKomendantovKononenko1996Neuron::runge_kutta_65 (void (*f) (double *, double *, double *, double), int dim, double dt, double * vars, double * params, double aux) {
    double apoyo[dim], retorno[dim];
    double k[6][dim];
    int j;

    (*f)(vars, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[0][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.2;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[1][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.075 + k[1][j] * 0.225;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[2][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.3 - k[1][j] * 0.9 + k[2][j] * 1.2;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[3][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.075 + k[1][j] * 0.675 - k[2][j] * 0.6 + k[3][j] * 0.75;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[4][j] = dt * retorno[j];
        apoyo[j] = vars[j]
                + k[0][j] * 0.660493827160493
                + k[1][j] * 2.5
                - k[2][j] * 5.185185185185185
                + k[3][j] * 3.888888888888889
                - k[4][j] * 0.864197530864197;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[5][j] = dt * retorno[j];
        apoyo[j] = vars[j]
                + k[0][j]*0.1049382716049382
                + k[2][j]*0.3703703703703703
                + k[3][j]*0.2777777777777777
                + k[4][j]*0.2469135802469135;
    }


    for(j = 0; j < dim; ++j) {
        vars[j] += k[0][j]*0.098765432098765+
                   k[2][j]*0.396825396825396+
                   k[3][j]*0.231481481481481+
                   k[4][j]*0.308641975308641-
                   k[5][j]*0.035714285714285;
    }

    return;
}

void RTHybridKomendantovKononenko1996Neuron::euler (void (*f) (double *, double *, double *, double), int dim, double dt, double * vars, double * params, double aux) {
    double apoyo[dim], retorno[dim];
    double k[4][dim];
    int j;

    (*f)(vars, retorno, params, aux);

    for(j = 0; j < dim; ++j) {
        vars[j] += dt * retorno[j];
    }

    return;
}


/**
 * @brief Selects the best integration step for a neuron model.
 * @param[in] dts Array with the integration steps
 * @param[in] pts Array with the numbers of points needed to complete one burst by the neuron model using the integration steps in dts
 * @param[in] length Length of the dts and pts arrays
 * @param[in] pts_live Number of points of burst from the living neuron
 * @param[out] dt Pointer to store the selected integration step-. If no integration step is found the value will be -1
 * @param[out] dt Pointer to store the number of points per burst from the neuron model with selected integration step. If no integration step is found the value will be -1
 */

void RTHybridKomendantovKononenko1996Neuron::select_dt_neuron_model (double * dts, double * pts, unsigned int length, double pts_live, double * dt, double * pts_burst) {
    double aux = pts_live;
    double factor = 1;
    double intpart, fractpart;
    int flag = 0;
    int i;

    *dt = -1;
    *pts_burst = -1;

    while (aux < pts[0]) {
        aux = pts_live * factor;
        factor += 1;

        for (i = length - 1; i >= 0; i--) {
        	if (pts[i] > aux) {
        		*dt = dts[i];
	            *pts_burst = pts[i];

	            fractpart = modf(*pts_burst / pts_live, &intpart);

	            if (fractpart <= 0.1*intpart) flag = 1;

	            break;
        	}
            
        }

        if (flag == 1) break;
    }

    if (flag == 0) {
        for (i = length - 1; i >= 0; i--) {
        	if (pts[i] > aux) {
        		*dt = dts[i];
	            *pts_burst = pts[i];

	            break;	
        	}
        }
    }

    return;
}


/**
 * @brief Sets Komendantov_Kononenko_1996 model number of points per burst and integration step.
 * 
 * If not previously specified by the user, the number of points per burst of the model and its integration step is set according to the living neuron number of points per burst.
 * @param[in] pts_match Number of points in a burst to match
 * @param[in] nm Pointer to the neuron model structure
 * @return Integration step, default -1
 */

double RTHybridKomendantovKononenko1996Neuron::set_pts_burst (double sec_per_burst) {
	int length = 0;
	int method = 3;
	double pts_match = sec_per_burst * freq;
	double pts_burst, dt;


	length = 19.000000;
	double dts[] = {0.000010, 0.000020, 0.000030, 0.000040, 0.000050, 0.000060, 0.000070, 0.000080, 0.000090, 0.000100, 0.000200, 0.000300, 0.000400, 0.000500, 0.000600, 0.000700, 0.000800, 0.000900, 0.001000};
	double pts[] = {489646.000000, 250715.000000, 172155.000000, 131029.714286, 105694.000000, 88541.272727, 76167.307692, 66740.857143, 59461.687500, 53611.000000, 26996.527778, 18041.654545, 13546.534247, 10844.543478, 9040.509091, 7750.382812, 6783.319728, 6030.793939, 5428.152174};

	select_dt_neuron_model(dts, pts, length, pts_match, &(params_model[NM_KOMENDANTOV_KONONENKO_1996_DT]), &(pts_burst));

	return pts_burst;
}


double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_V (double * vars, double * params) {
  double i_na_ttx = nm_komendantov_kononenko_1996_i_na_ttx(vars, params);
  double i_k_tea = nm_komendantov_kononenko_1996_i_k_tea(vars, params);
  double i_k = nm_komendantov_kononenko_1996_i_k(vars, params);
  double i_na = nm_komendantov_kononenko_1996_i_na(vars, params);
  double i_na_v = nm_komendantov_kononenko_1996_i_na_v(vars, params);
  double i_b = nm_komendantov_kononenko_1996_i_b(vars, params);
  double i_ca = nm_komendantov_kononenko_1996_i_ca(vars, params);
  double i_ca_ca = nm_komendantov_kononenko_1996_i_ca_ca(vars, params);

  return (-(i_na_ttx + i_k_tea + i_k + i_na + i_na_v + i_b + i_ca + i_ca_ca) + params[NM_KOMENDANTOV_KONONENKO_1996_I] - params[NM_KOMENDANTOV_KONONENKO_1996_SYN]) / params[NM_KOMENDANTOV_KONONENKO_1996_CM];
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_i_na_v (double * vars, double * params) {
  return params[NM_KOMENDANTOV_KONONENKO_1996_G_NA_V] * (1.000000 / (1.000000 + exp(-0.200000 * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] + 45.000000)))) * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] - params[NM_KOMENDANTOV_KONONENKO_1996_V_NA]);
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_i_k (double * vars, double * params) {
  return params[NM_KOMENDANTOV_KONONENKO_1996_G_K] * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] - params[NM_KOMENDANTOV_KONONENKO_1996_V_K]);
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_i_na (double * vars, double * params) {
  return params[NM_KOMENDANTOV_KONONENKO_1996_G_NA] * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] - params[NM_KOMENDANTOV_KONONENKO_1996_V_NA]);
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_i_b (double * vars, double * params) {
  return params[NM_KOMENDANTOV_KONONENKO_1996_G_B] * vars[NM_KOMENDANTOV_KONONENKO_1996_M_B] * vars[NM_KOMENDANTOV_KONONENKO_1996_H_B] * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] - params[NM_KOMENDANTOV_KONONENKO_1996_V_B]);
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_m_b (double * vars, double * params) {
  return (1.000000 / (1.000000 + exp(0.400000 * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] + 34.000000))) - vars[NM_KOMENDANTOV_KONONENKO_1996_M_B]) / 0.050000;
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_h_b (double * vars, double * params) {
  return (1.000000 / (1.000000 + exp(-0.550000 * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] + 43.000000))) - vars[NM_KOMENDANTOV_KONONENKO_1996_H_B]) / 1.500000;
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_i_na_ttx (double * vars, double * params) {
  return params[NM_KOMENDANTOV_KONONENKO_1996_G_NA_TTX] * vars[NM_KOMENDANTOV_KONONENKO_1996_M]*vars[NM_KOMENDANTOV_KONONENKO_1996_M]*vars[NM_KOMENDANTOV_KONONENKO_1996_M] * vars[NM_KOMENDANTOV_KONONENKO_1996_H] * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] - params[NM_KOMENDANTOV_KONONENKO_1996_V_NA]);
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_i_k_tea (double * vars, double * params) {
  return params[NM_KOMENDANTOV_KONONENKO_1996_G_K_TEA] * vars[NM_KOMENDANTOV_KONONENKO_1996_N]*vars[NM_KOMENDANTOV_KONONENKO_1996_N]*vars[NM_KOMENDANTOV_KONONENKO_1996_N]*vars[NM_KOMENDANTOV_KONONENKO_1996_N] * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] - params[NM_KOMENDANTOV_KONONENKO_1996_V_K]);
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_m (double * vars, double * params) {
  //printf("(1 / 1 + exp(-0.4 * (V[%f] + 31))) - m[%f]) / 0.0005\n", vars[NM_KOMENDANTOV_KONONENKO_1996_V], vars[NM_KOMENDANTOV_KONONENKO_1996_M]);
  return (1.0 / (1.0 + exp(-0.4 * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] + 31.0))) - vars[NM_KOMENDANTOV_KONONENKO_1996_M]) / 0.0005;
  //return (1.000000 / (1.000000 + exp(-0.400000 * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] + 31.000000))) - vars[NM_KOMENDANTOV_KONONENKO_1996_M]) / 0.000500;
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_h (double * vars, double * params) {
  return (1.000000 / (1.000000 + exp(0.250000 * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] + 45.000000))) - vars[NM_KOMENDANTOV_KONONENKO_1996_H]) / 0.010000;
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_n (double * vars, double * params) {
  return (1.000000 / (1.000000 + exp(-0.180000 * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] + 25.000000))) - vars[NM_KOMENDANTOV_KONONENKO_1996_N]) / 0.015000;
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_i_ca (double * vars, double * params) {
  return params[NM_KOMENDANTOV_KONONENKO_1996_G_CA] * vars[NM_KOMENDANTOV_KONONENKO_1996_M_CA]*vars[NM_KOMENDANTOV_KONONENKO_1996_M_CA] * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] - params[NM_KOMENDANTOV_KONONENKO_1996_V_CA]);
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_m_ca (double * vars, double * params) {
  return (1.000000 / (1.000000 + exp(-0.200000 * vars[NM_KOMENDANTOV_KONONENKO_1996_V])) - vars[NM_KOMENDANTOV_KONONENKO_1996_M_CA]) / 0.010000;
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_i_ca_ca (double * vars, double * params) {
  return params[NM_KOMENDANTOV_KONONENKO_1996_G_CA_CA] * (1.000000 / (1.000000 + exp(-0.060000 * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] + 45.000000)))) * (1.000000 / (1.000000 + exp(params[NM_KOMENDANTOV_KONONENKO_1996_K_BETA] * (vars[NM_KOMENDANTOV_KONONENKO_1996_CA] - params[NM_KOMENDANTOV_KONONENKO_1996_BETA])))) * (vars[NM_KOMENDANTOV_KONONENKO_1996_V] - params[NM_KOMENDANTOV_KONONENKO_1996_V_CA]);
}

double RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_Ca (double * vars, double * params) {
  double i_ca = nm_komendantov_kononenko_1996_i_ca(vars, params);
  return params[NM_KOMENDANTOV_KONONENKO_1996_RHO] * ((-i_ca / 808.310846) - (params[NM_KOMENDANTOV_KONONENKO_1996_K_S] * vars[NM_KOMENDANTOV_KONONENKO_1996_CA]));
}

void RTHybridKomendantovKononenko1996Neuron::nm_komendantov_kononenko_1996_f (double * vars, double * ret, double * params, double syn) {
  params[NM_KOMENDANTOV_KONONENKO_1996_SYN] = syn;

  ret[NM_KOMENDANTOV_KONONENKO_1996_V] = nm_komendantov_kononenko_1996_V(vars, params);
  ret[NM_KOMENDANTOV_KONONENKO_1996_M_B] = nm_komendantov_kononenko_1996_m_b(vars, params);
  ret[NM_KOMENDANTOV_KONONENKO_1996_H_B] = nm_komendantov_kononenko_1996_h_b(vars, params);
  ret[NM_KOMENDANTOV_KONONENKO_1996_M] = nm_komendantov_kononenko_1996_m(vars, params);
  ret[NM_KOMENDANTOV_KONONENKO_1996_H] = nm_komendantov_kononenko_1996_h(vars, params);
  ret[NM_KOMENDANTOV_KONONENKO_1996_N] = nm_komendantov_kononenko_1996_n(vars, params);
  ret[NM_KOMENDANTOV_KONONENKO_1996_M_CA] = nm_komendantov_kononenko_1996_m_ca(vars, params);
  ret[NM_KOMENDANTOV_KONONENKO_1996_CA] = nm_komendantov_kononenko_1996_Ca(vars, params);
}


void
RTHybridKomendantovKononenko1996Neuron::execute(void)
{
  int i;

  if (burst_duration_value <= -1) {
      burst_duration = input(1);
      s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq));
      if (s_points < 1) s_points = 1;
  }
  
  for (i = 0; i < s_points; i++) {
  	runge_kutta_65(nm_komendantov_kononenko_1996_f, NM_KOMENDANTOV_KONONENKO_1996_N_VARS, params_model[NM_KOMENDANTOV_KONONENKO_1996_DT], vars_model, params_model, input(0));	
  }


  output(0) = vars_model[NM_KOMENDANTOV_KONONENKO_1996_V] / 1000.0;
  output(1) = vars_model[NM_KOMENDANTOV_KONONENKO_1996_V];
  return;
}

void
RTHybridKomendantovKononenko1996Neuron::initParameters(void)
{
  burst_duration_value = 1.0;
  burst_duration = burst_duration_value;
  freq = 1.0 / (period * 1e-3);
  s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq)); 
  if (s_points == 0) s_points = 1;


  params_model[NM_KOMENDANTOV_KONONENKO_1996_I] = 0;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_SYN] = 0;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_CM] = 0.02;

  params_model[NM_KOMENDANTOV_KONONENKO_1996_G_NA_V] = 0.11;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_V_NA] = 40;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_G_NA] = 0.0231;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_G_NA_TTX] = 400;

  params_model[NM_KOMENDANTOV_KONONENKO_1996_G_K] = 0.25;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_G_K_TEA] = 10;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_V_K] = -70;

  params_model[NM_KOMENDANTOV_KONONENKO_1996_G_B] = 0.165;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_V_B] = -58;

  params_model[NM_KOMENDANTOV_KONONENKO_1996_G_CA] = 1.5;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_V_CA] = 150;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_G_CA_CA] = 0.02;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_K_BETA] = 15000;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_BETA] = 0.00004;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_RHO] = 0.002;
  params_model[NM_KOMENDANTOV_KONONENKO_1996_K_S] = 50;

  vars_model[NM_KOMENDANTOV_KONONENKO_1996_V] = -55;
  vars_model[NM_KOMENDANTOV_KONONENKO_1996_CA] = 0;
  vars_model[NM_KOMENDANTOV_KONONENKO_1996_M_B] = (1.000000 / (1.000000 + exp(0.400000 * (vars_model[NM_KOMENDANTOV_KONONENKO_1996_V] + 34.000000))));
  vars_model[NM_KOMENDANTOV_KONONENKO_1996_H_B] = (1.000000 / (1.000000 + exp(-0.550000 * (vars_model[NM_KOMENDANTOV_KONONENKO_1996_V] + 43.000000))));
  vars_model[NM_KOMENDANTOV_KONONENKO_1996_M] = (1.000000 / (1.000000 + exp(-0.400000 * (vars_model[NM_KOMENDANTOV_KONONENKO_1996_V] + 31.000000))));
  vars_model[NM_KOMENDANTOV_KONONENKO_1996_H] = (1.000000 / (1.000000 + exp(0.250000 * (vars_model[NM_KOMENDANTOV_KONONENKO_1996_V] + 45.000000))));
  vars_model[NM_KOMENDANTOV_KONONENKO_1996_N] = (1.000000 / (1.000000 + exp(-0.180000 * (vars_model[NM_KOMENDANTOV_KONONENKO_1996_V] + 25.000000))));
  vars_model[NM_KOMENDANTOV_KONONENKO_1996_M_CA] = (1.000000 / (1.000000 + exp(-0.200000 * vars_model[NM_KOMENDANTOV_KONONENKO_1996_V])));

}

void
RTHybridKomendantovKononenko1996Neuron::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag) {
    case INIT:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      freq = 1.0 / (period * 1e-3);

      setParameter("Burst duration (s)", burst_duration_value);
      setParameter("i", params_model[NM_KOMENDANTOV_KONONENKO_1996_I]);
      setParameter("cm", params_model[NM_KOMENDANTOV_KONONENKO_1996_CM]);

      setParameter("g_na_v", params_model[NM_KOMENDANTOV_KONONENKO_1996_G_NA_V]);
      setParameter("v_na", params_model[NM_KOMENDANTOV_KONONENKO_1996_V_NA]);
      setParameter("g_na", params_model[NM_KOMENDANTOV_KONONENKO_1996_G_NA]);
      setParameter("g_na_ttx", params_model[NM_KOMENDANTOV_KONONENKO_1996_G_NA_TTX]);

      setParameter("g_k", params_model[NM_KOMENDANTOV_KONONENKO_1996_G_K]);
      setParameter("g_k_tea", params_model[NM_KOMENDANTOV_KONONENKO_1996_G_K_TEA]);
      setParameter("v_k", params_model[NM_KOMENDANTOV_KONONENKO_1996_V_K]);

      setParameter("g_b", params_model[NM_KOMENDANTOV_KONONENKO_1996_G_B]);
      setParameter("v_b", params_model[NM_KOMENDANTOV_KONONENKO_1996_V_B]);

      setParameter("g_ca", params_model[NM_KOMENDANTOV_KONONENKO_1996_G_CA]);
      setParameter("v_ca", params_model[NM_KOMENDANTOV_KONONENKO_1996_V_CA]);
      setParameter("g_ca_ca", params_model[NM_KOMENDANTOV_KONONENKO_1996_G_CA_CA]);
      setParameter("k_beta", params_model[NM_KOMENDANTOV_KONONENKO_1996_K_BETA]);
      setParameter("beta", params_model[NM_KOMENDANTOV_KONONENKO_1996_BETA]);
      setParameter("rho", params_model[NM_KOMENDANTOV_KONONENKO_1996_RHO]);
      setParameter("k_s", params_model[NM_KOMENDANTOV_KONONENKO_1996_K_S]);

      setParameter("v0 (mV)", vars_model[NM_KOMENDANTOV_KONONENKO_1996_V]);

      setState("v", vars_model[NM_KOMENDANTOV_KONONENKO_1996_V]);
      setState("s_points", s_points);
      setState("dt", params_model[NM_KOMENDANTOV_KONONENKO_1996_DT]);
      setState("syn", params_model[NM_KOMENDANTOV_KONONENKO_1996_SYN]);

      break;

    case MODIFY:
      burst_duration_value = getParameter("Burst duration (s)").toDouble();
      freq = 1.0 / (period * 1e-3);

      if (burst_duration_value <= -1) {
          burst_duration = input(1);
      } else {
          burst_duration = burst_duration_value;
      }
      
      s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq));
      if (s_points < 1) s_points = 1;

      params_model[NM_KOMENDANTOV_KONONENKO_1996_I] = getParameter("i").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_CM] = getParameter("cm").toDouble();

      params_model[NM_KOMENDANTOV_KONONENKO_1996_G_NA_V] = getParameter("g_na_v").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_V_NA] = getParameter("v_na").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_G_NA] = getParameter("g_na").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_G_NA_TTX] = getParameter("g_na_ttx").toDouble();

      params_model[NM_KOMENDANTOV_KONONENKO_1996_G_K] = getParameter("g_k").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_G_K_TEA] = getParameter("g_k_tea").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_V_K] = getParameter("v_k").toDouble();

      params_model[NM_KOMENDANTOV_KONONENKO_1996_G_B] = getParameter("g_b").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_V_B] = getParameter("v_b").toDouble();

      params_model[NM_KOMENDANTOV_KONONENKO_1996_G_CA] = getParameter("g_ca").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_V_CA] = getParameter("v_ca").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_G_CA_CA] = getParameter("g_ca_ca").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_K_BETA] = getParameter("k_beta").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_BETA] = getParameter("beta").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_RHO] = getParameter("rho").toDouble();
      params_model[NM_KOMENDANTOV_KONONENKO_1996_K_S] = getParameter("k_s").toDouble();

      vars_model[NM_KOMENDANTOV_KONONENKO_1996_V] = getParameter("v0 (mV)").toDouble();
      break;

    case UNPAUSE:
      freq = 1.0 / (period * 1e-3);
      s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq)); 
      if (s_points == 0) s_points = 1;

      break;

    case PAUSE:
      output(0) = 0;
      output(1) = 0; 
      break;

    case PERIOD:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      freq = 1.0 / (period * 1e-3);
      s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq)); 
      if (s_points == 0) s_points = 1;
      break;

    default:
      break;
  }
}

void
RTHybridKomendantovKononenko1996Neuron::customizeGUI(void)
{
  QGridLayout* customlayout = DefaultGUIModel::getLayout();

  QGroupBox* button_group = new QGroupBox;

  QPushButton* abutton = new QPushButton("Button A");
  QPushButton* bbutton = new QPushButton("Button B");
  QHBoxLayout* button_layout = new QHBoxLayout;
  button_group->setLayout(button_layout);
  button_layout->addWidget(abutton);
  button_layout->addWidget(bbutton);
  QObject::connect(abutton, SIGNAL(clicked()), this, SLOT(aBttn_event()));
  QObject::connect(bbutton, SIGNAL(clicked()), this, SLOT(bBttn_event()));

  customlayout->addWidget(button_group, 0, 0);
  setLayout(customlayout);
}

// functions designated as Qt slots are implemented as regular C++ functions
void
RTHybridKomendantovKononenko1996Neuron::aBttn_event(void)
{
}

void
RTHybridKomendantovKononenko1996Neuron::bBttn_event(void)
{
}
