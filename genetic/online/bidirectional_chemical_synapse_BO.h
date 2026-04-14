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

#ifndef RTHYBRID_BIDIRECTIONAL_CHEMICAL_SYNAPSE_BO_H
#define RTHYBRID_BIDIRECTIONAL_CHEMICAL_SYNAPSE_BO_H

#include <default_gui_model.h>
#include <kfr/all.hpp>
#include <thread>
#include <atomic>
#include "utils.hpp"

// RTXI module implementing a bidirectional chemical synapse model.
// A non-real-time (NRT) thread runs Bayesian Optimization (BO) and coordinates
// with the real-time (RT) `execute()` loop via lock-free atomics and
// double-buffered parameter slots.
class BidirectionalChemicalSynapseBO : public DefaultGUIModel
{

  Q_OBJECT

public:
  // Standard RTXI plugin lifecycle.
  BidirectionalChemicalSynapseBO(void);
  virtual ~BidirectionalChemicalSynapseBO(void);

  // RTXI real-time entry point (runs once per RT cycle).
  void execute(void);

  // DefaultGUIModel GUI helpers.
  void createGUI(DefaultGUIModel::variable_t *, int);
  void customizeGUI(void);

  // Limbo evaluation functor needs access to private helpers/state.
  friend struct EvaluationFunctor;

protected:
  virtual void update(DefaultGUIModel::update_flags_t);

private:
  // BO runtime configuration (GUI parameters/state).
  // - evaluations_completed: number of candidate evaluations completed (STATE)
  // - initial_samples: LHS initialization samples before BO iterations
  // - iterations: BO iterations after initialization
  // - evaluation_time: signal recording duration per candidate (ms)
  // - stabilization_time: wait after swapping params before recording (ms)
  // - search_phase: affects E_syn range and shape-score polarity (1/0)
  double evaluations_completed;
  unsigned int initial_samples;
  unsigned int iterations;
  double evaluation_time, stabilization_time;
  unsigned int search_phase;

  // Target current range to achieve during BO (per direction, in nA).
  double expected_i_min_12, expected_i_max_12, expected_i_min_21, expected_i_max_21;

  // Cutoff frequencies (kHz) used to separate I_fast vs I_slow for scoring.
  double fc_1, fc_2;

  // Voltage bounds for parameter-range initialization (can be fixed or dynamic).
  unsigned int dynamic_v_min_max_1, dynamic_v_min_max_2;
  double v_min_1, v_max_1, v_min_2, v_max_2;

  // Output clamps for total current (safety bounds), per direction (nA).
  double i_min_12, i_max_12, i_min_21, i_max_21;

  // Integration timestep used by the slow synaptic gating ODE.
  double dt, dt_factor;

  // Select whether to include fast and/or slow current components in the model.
  unsigned int use_i_fast_12, use_i_slow_12, use_i_fast_21, use_i_slow_21;

  // Double-buffered synapse parameters: RT reads one slot while NRT writes the other.
  ChemicalSynapseParams params_21[2];
  ChemicalSynapseParams params_12[2];

  // Slow gating-variable state (maintained in RT loop).
  double m_slow_21, m_slow_12;
  double period;

  // BO thread and stop flag (set from GUI/parameter changes).
  std::thread BO_NRT_thread;
  std::atomic<bool> stop_BO;
  bool BO_running;

  // Synchronization between RT and NRT:
  // - synapse_idx: which params slot RT should read (release written by NRT)
  // - last_synapse_idx_read_RT: which slot RT last observed (handshake for NRT)
  std::atomic<size_t> synapse_idx;
  std::atomic<size_t> last_synapse_idx_read_RT;

  // RT->NRT signal capture handshake.
  // - RT_storing: set by NRT to request capture; cleared by RT when done.
  // - storing_idx/num_elements: capture progress and requested length.
  std::atomic<bool> RT_storing;
  size_t storing_idx;
  size_t num_elements;

  // Signal buffers filled in RT and scored in NRT (sized during BO start).
  kfr::univector<double> i_fast_sig_12;
  kfr::univector<double> i_fast_sig_21;
  kfr::univector<double> i_slow_sig_12;
  kfr::univector<double> i_slow_sig_21;
  kfr::univector<double> v_sig_1;
  kfr::univector<double> v_sig_2;

  // GUI control for starting/stopping BO.
  QPushButton *BO_button;

  // Initialize defaults and internal state.
  void initParameters();

  // 6-stage Runge–Kutta integrator used for the slow gating ODE (fixed step dt).
  void runge_kutta_65(double (*f)(double, double, const ChemicalSynapseParams &), double &m_slow, double v_pre, double dt, const ChemicalSynapseParams &params);

  // Compute synaptic currents (slow component updates m_slow; fast is instantaneous).
  double compute_i_slow(double &m_slow, double v_pre, double v_post, const ChemicalSynapseParams &params);
  double compute_i_fast(double v_pre, double v_post, const ChemicalSynapseParams &params);

  // ODE RHS for the slow gating variable.
  static double sm_chemical_synapse_m(double m_slow, double v_pre, const ChemicalSynapseParams &params);

  // Evaluate a BO candidate by swapping params into RT, capturing signals, and scoring.
  ChemicalSynapseEvaluation evaluate_candidate(
      const Candidate &candidate,
      double fs,
      size_t effective_pad_12,
      size_t effective_pad_21,
      EvaluationPadBuffers &pad_buffers,
      double max_i_dist_12,
      double max_i_dist_21,
      size_t &curr_synapse_idx);

  // Decode normalized BO vector x\in[0,1]^d into physical synapse parameters.
  Candidate decode_to_candidate(const Eigen::VectorXd &x,
                                const BOParamRanges &ranges_12,
                                const BOParamRanges &ranges_21);

  // Non-real-time BO loop (runs in its own thread).
  void NRT_BO(double period_t);

  // Lock/unlock GUI fields while BO is running.
  void set_params_read_only(bool read_only);

  // Helpers for initializing parameter structs and coordinating RT/NRT handshakes.
  void init_syn_params_and_vars(ChemicalSynapseParams &params);
  bool wait_until_RT_read_idx_or_stop(size_t idx_to_achieve);

private slots:
  void toggle_BO_event(void);
  void stop_BO_event_async(void);
  void update_params_gui(void);
  void set_evaluations_completed(double evals);
};

#endif