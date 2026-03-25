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

class BinarySemaphore
{
private:
  std::mutex mtx_;
  std::condition_variable cv_;
  bool available_;

public:
  explicit BinarySemaphore(bool initial_state = false) : available_(initial_state) {}

  void acquire()
  {
    std::unique_lock<std::mutex> lock(mtx_);
    cv_.wait(lock, [this]()
             { return available_; });
    available_ = false;
  }

  bool try_acquire()
  {
    std::lock_guard<std::mutex> lock(mtx_);
    if (available_)
    {
      available_ = false;
      return true;
    }
    return false;
  }

  void release()
  {
    {
      std::lock_guard<std::mutex> lock(mtx_);
      available_ = true;
    }
    cv_.notify_one();
  }
};

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
  double burst_duration, burst_duration_gui;
  double scale_21, offset_21, scale_12, offset_12;
  double scale_21_gui, scale_12_gui;
  unsigned int is_living_1, is_living_2;
  unsigned int dynamic_offset_21, dynamic_offset_12;
  double s_points;

  double evaluation_time;
  double stabilization_time;
  unsigned int search_phase;
  double current_max;
  double current_min;
  std::thread genetic_NRT_thread;
  std::atomic<bool> RT_storing;
  std::atomic<bool> stop_genetic;
  std::atomic<bool> genetic_running;
  BinarySemaphore synapse_lock{true};
  BinarySemaphore scaling_factors_21_lock{true}, scaling_factors_12_lock{true};
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

  void select_dt_neuron_model(double *dts, double *pts, size_t length, double pts_live, double *dt, double *pts_burst);
  double set_pts_burst(double sec_per_burst);

  static double sm_chemical_synapse_m(double m_slow, double v_pre, const ChemicalSynapseParams &params);

  void NRT_genetic(double thread_freq, double scale_21_thread, double offset_21_thread, double scale_12_thread, double offset_12_thread);

  void set_params_read_only(bool read_only);
  void set_param_read_only(const QString &name, const QPalette &pal, bool read_only);

private slots:
  void toggle_genetic_event(void);
  void stop_genetic_event_async(void);
  void update_params_gui(void);
};

#endif