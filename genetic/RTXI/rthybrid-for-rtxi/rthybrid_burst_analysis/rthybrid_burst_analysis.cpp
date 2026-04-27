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

#include "rthybrid_burst_analysis.h"
#include <iostream>
#include <main_window.h>

extern "C" Plugin::Object*
createRTXIPlugin(void)
{
  return new RTHybridBurstAnalysis();
}

static DefaultGUIModel::variable_t vars[] = {
    { "Observation time (s)", "Duration of the interval where the minimum and maximum membrane potential and the burst duration will be detected (in s)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

    { "Vm (V)", "Membrane potential (in V)", DefaultGUIModel::INPUT,},

    { "Min (V)", "Minimum membrane potential (in V)", DefaultGUIModel::OUTPUT,},
    { "Max (V)", "Maximum membrane potential (in V)", DefaultGUIModel::OUTPUT,},
    { "Burst duration (s)", "Average burst duration (in s)", DefaultGUIModel::OUTPUT, },

    { "Min (V)", "Minimum membrane potential (in V)", DefaultGUIModel::STATE,},
    { "Max (V)", "Maximum membrane potential (in V)", DefaultGUIModel::STATE,},
    { "Burst duration (s)", "Average burst duration (in s)", DefaultGUIModel::STATE, },
    { "Upper threshold (V)", "Upper threshold for burst detection (in V)", DefaultGUIModel::STATE,},
    { "Lower threshold (V)", "Lower threshold for burst detection (in V)", DefaultGUIModel::STATE,},
    { "Amplitude (V)", "Membrane potential amplitude (in V)", DefaultGUIModel::STATE,},
    { "Points counter", "Number of points so far in the current interval", DefaultGUIModel::STATE,},
    { "Burst counter", "Number of bursts so far in the current interval", DefaultGUIModel::STATE,},
    { "Is burst", "1 if it is currently a burst or 0 if not", DefaultGUIModel::STATE,},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

RTHybridBurstAnalysis::RTHybridBurstAnalysis(void)
  : DefaultGUIModel("RTHybridBurstAnalysis", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>RTHybridBurstAnalysis:</b><br>RTHybrid module for RTXI to get the minimum and maximum membrane potential values of a neuron and its bursts duration.</p>");
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

RTHybridBurstAnalysis::~RTHybridBurstAnalysis(void)
{
}

void
RTHybridBurstAnalysis::execute(void)
{
    double vm = input(0); //V
    double time = count * period; //ms

    if (time > (observation_time * 1000)) {
        min = temp_min;
        max = temp_max;
        output(0) = min;
        output(1) = max;
        temp_min = 999999;
        temp_max = -999999;
        count = 0;

        sec_per_burst = (burst_dur_sum / burst_counter) / freq;
        output(2) = sec_per_burst;
        pts_counter = 0;
        burst_counter = 0;
        burst_dur_sum = 0;
        old_burst_time = 0;
    }

    if (vm < temp_min) temp_min = vm;
    if (vm > temp_max) temp_max = vm;

    range = max - min;
    thresh_down = min + (range * 0.1);
    thresh_up = min + (range * 0.9);


    if (is_burst == 0 && vm > thresh_up) {
        is_burst = 1;
        burst_counter++;
        burst_dur_sum += pts_counter - old_burst_time;
        old_burst_time = pts_counter;
    } else if (is_burst == 1 && vm < thresh_down) {
        is_burst = 0;
    }

    pts_counter++;
    count++;
}

void
RTHybridBurstAnalysis::initParameters(void)
{
    observation_time = 5;

    min = 999999;
    max = -999999;
    temp_min = 999999;
    temp_max = -999999;
    count = 0;

    sec_per_burst = 0;
    pts_counter = 0;
    burst_counter = -1;
    is_burst = 1;
    burst_dur_sum = 0;
    old_burst_time = 0;
}

void
RTHybridBurstAnalysis::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag) {
    case INIT:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      freq = 1.0 / (period * 1e-3);
      setParameter("Observation time (s)", observation_time);

      setState("Min (V)", min);
      setState("Max (V)", max);
      setState("Burst duration (s)", sec_per_burst);

      setState("Upper threshold (V)", thresh_up);
      setState("Lower threshold (V)", thresh_down);
      setState("Amplitude (V)", range);
      setState("Is burst", is_burst);
      setState("Points counter", pts_counter);
      setState("Burst counter", burst_counter);

      break;

    case MODIFY:
      observation_time = getParameter("Observation time (s)").toDouble(); //s
      break;

    case UNPAUSE:
      break;

    case PAUSE:
      break;

    case PERIOD:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      freq = 1.0 / (period * 1e-3);
      break;

    default:
      break;
  }
}

void
RTHybridBurstAnalysis::customizeGUI(void)
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
RTHybridBurstAnalysis::aBttn_event(void)
{
}

void
RTHybridBurstAnalysis::bBttn_event(void)
{
}
