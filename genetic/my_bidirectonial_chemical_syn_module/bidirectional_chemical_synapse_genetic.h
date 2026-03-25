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

#ifndef RTHYBRID_CHEMICAL_SYNAPSE_GENETIC_H
#define RTHYBRID_CHEMICAL_SYNAPSE_GENETIC_H

#include <default_gui_model.h>

struct ChemicalSynapseParams
{
  double e_syn;
  double g_fast;
  double s_fast;
  double v_fast;
  double g_slow;
  double k1;
  double k2;
  double s_slow;
  double v_slow;
};

class BidirectionalChemicalSynapseGenetic : public DefaultGUIModel
{

  Q_OBJECT

public:
  BidirectionalChemicalSynapseGenetic(void);
  virtual ~BidirectionalChemicalSynapseGenetic(void);

  void execute(void);
  void createGUI(DefaultGUIModel::variable_t *, int);
  void customizeGUI(void);

protected:
  virtual void update(DefaultGUIModel::update_flags_t);

private:
  double m_slow_21, m_slow_12;
  ChemicalSynapseParams params_21;
  ChemicalSynapseParams params_12;
  unsigned int use_i_fast_21, use_i_slow_21, use_i_fast_12, use_i_slow_12;
  double dt;
  double period, freq;
  double last_burst_duration, burst_duration_gui;
  double offset_21, offset_12;
  double scale_21_gui, scale_12_gui;
  unsigned int dynamic_offset_21, dynamic_offset_12;
  double s_points;

  void initParameters();

  void runge_kutta_65(double (*f)(double, double, const ChemicalSynapseParams &), double &m_slow, double v_pre, double dt, const ChemicalSynapseParams &params);
  double compute_synapse_current(double &m_slow, double v_pre, double v_post, const ChemicalSynapseParams &params, unsigned int use_i_fast, unsigned int use_i_slow);
  void select_dt_neuron_model(double *dts, double *pts, size_t length, double pts_live, double *dt, double *pts_burst);
  double set_pts_burst(double sec_per_burst);

  static double sm_chemical_synapse_m(double m_slow, double v_pre, const ChemicalSynapseParams &params);

private slots:
  void aBttn_event(void);
  void bBttn_event(void);
};

#endif