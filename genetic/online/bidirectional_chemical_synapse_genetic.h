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
  double m_slow_21[2], m_slow_12[2];
  ChemicalSynapseParams params_21[2];
  ChemicalSynapseParams params_12[2];
  unsigned int use_i_fast_21, use_i_slow_21, use_i_fast_12, use_i_slow_12;
  double dt;
  double period, freq;
  double burst_duration, burst_duration_gui;
  double scale_21[2], offset_21[2], scale_12[2], offset_12[2];
  double scale_21_gui, scale_12_gui;
  unsigned int dynamic_offsets;
  unsigned int dynamic_min_max_1;
  double max_1, min_1;
  double s_points;

  double evaluation_time;
  double stabilization_time;
  unsigned int num_generations, generations_completed;
  unsigned int population_size, individuals_completed;
  unsigned int search_phase;
  unsigned int i_ranges_from_neuron;
  double i_max_21;
  double i_min_21;
  double i_max_12;
  double i_min_12;
  std::thread genetic_NRT_thread;
  std::atomic<bool> RT_storing;
  std::atomic<bool> stop_genetic;
  bool genetic_running;
  std::atomic<size_t> synapse_idx;
  std::atomic<size_t> scaling_factors_21_idx;
  std::atomic<size_t> scaling_factors_12_idx;
  size_t storing_idx;
  size_t num_elements;

  univector<double> i_fast_sig_12;
  univector<double> i_fast_sig_21;
  univector<double> i_slow_12;
  univector<double> i_slow_21;
  univector<double> v1_scaled_sig;
  univector<double> v1_sig;
  univector<double> v2_scaled_sig;
  univector<double> v2_sig;

  QPushButton *gentic_button;

  void initParameters();

  void runge_kutta_65(double (*f)(double, double, const ChemicalSynapseParams &), double &m_slow, double v_pre, double dt, const ChemicalSynapseParams &params);

  double compute_i_slow(double &m_slow, double v_pre, double v_post, const ChemicalSynapseParams &params);
  double compute_i_fast(double v_pre, double v_post, const ChemicalSynapseParams &params);

  void select_dt_neuron_model(const double *dts, const double *pts, size_t length, double pts_live, double *dt, double *pts_burst);
  double set_pts_burst(double sec_per_burst);

  static double sm_chemical_synapse_m(double m_slow, double v_pre, const ChemicalSynapseParams &params);

  void NRT_genetic(double thread_freq, double scale_21_t, double offset_21_t, double scale_12_t, double offset_12_t, unsigned int i_ranges_from_neuron_t, double i_max_21_t, double i_min_21_t, double i_max_12_t, double i_min_12_t, double max_1_t, double min_1_t);

  void set_params_read_only(bool read_only);

private slots:
  void toggle_genetic_event(void);
  void stop_genetic_event_async(void);
  void update_params_gui(void);
  void set_generations_completed(size_t generations);
  void set_individuals_completed(size_t individuals);
};

#endif