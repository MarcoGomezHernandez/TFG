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
  double vars_model[2];
  double params_model[25];
  double period, freq;
  double burst_duration, burst_duration_value;
  double scale21, scale21_value;
  double offset21, offset21_value;
  double scale12, scale12_value;
  double offset12, offset12_value;
  double s_points;

  void initParameters();

  void runge_kutta_65(void (*f)(double *, double *, double *), int dim, double dt, double *vars, double *params);
  void select_dt_neuron_model(double *dts, double *pts, unsigned int length, double pts_live, double *dt, double *pts_burst);
  double set_pts_burst(double sec_per_burst);

  static void sm_chemical_synapse_m_21(double *vars, double *ret, double *params);
  static void sm_chemical_synapse_m_12(double *vars, double *ret, double *params);

private slots:
  void aBttn_event(void);
  void bBttn_event(void);
};

#endif