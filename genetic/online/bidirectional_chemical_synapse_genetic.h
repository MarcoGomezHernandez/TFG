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
#include <kfr/all.hpp>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>

using namespace kfr;

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
  double generations_completed, individuals_completed;
  unsigned int num_generations, population_size;
  double evaluation_time, stabilization_time;
  unsigned int search_phase;
  double i_max_12, i_min_12, i_max_21, i_min_21;

  unsigned int dynamic_v_min_max_1, dynamic_v_min_max_2;
  double v_max_1, v_min_1, v_max_2, v_min_2;

  double dt, dt_factor;

  unsigned int use_i_fast_12, use_i_slow_12, use_i_fast_21, use_i_slow_21;

  ChemicalSynapseParams params_21[2];
  ChemicalSynapseParams params_12[2];

  double m_slow_21, m_slow_12;
  double period;

  std::thread genetic_NRT_thread;
  std::atomic<bool> stop_genetic;
  bool genetic_running;

  std::atomic<size_t> synapse_idx;

  std::atomic<bool> RT_storing;
  size_t storing_idx;
  size_t num_elements;

  univector<double> i_fast_sig_12;
  univector<double> i_fast_sig_21;
  univector<double> i_slow_sig_12;
  univector<double> i_slow_sig_21;
  univector<double> v_sig_1;
  univector<double> v_sig_2;

  QPushButton *gentic_button;

  void initParameters();

  void runge_kutta_65(double (*f)(double, double, const ChemicalSynapseParams &), double &m_slow, double v_pre, double dt, const ChemicalSynapseParams &params);

  double compute_i_slow(double &m_slow, double v_pre, double v_post, const ChemicalSynapseParams &params);
  double compute_i_fast(double v_pre, double v_post, const ChemicalSynapseParams &params);

  static double sm_chemical_synapse_m(double m_slow, double v_pre, const ChemicalSynapseParams &params);

  void NRT_genetic(double period_t, double dt_factor_t, double v_max_1_t, double v_min_1_t, double v_max_2_t, double v_min_2_t);

  void set_params_read_only(bool read_only);

  void init_syn_params_and_vars(ChemicalSynapseParams &params);

private slots:
  void toggle_genetic_event(void);
  void stop_genetic_event_async(void);
  void update_params_gui(void);
  void set_generations_completed(double generations);
  void set_individuals_completed(double individuals);
};

#endif