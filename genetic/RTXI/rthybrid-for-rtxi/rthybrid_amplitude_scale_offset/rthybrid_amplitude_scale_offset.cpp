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

#include "rthybrid_amplitude_scale_offset.h"
#include <iostream>
#include <main_window.h>

extern "C" Plugin::Object*
createRTXIPlugin(void)
{
  return new RTHybridAmplitudeScaleOffset();
}

static DefaultGUIModel::variable_t vars[] = {
  {"Min 1 (V)", "Neuron 1 minimum membrane potential (in V)", DefaultGUIModel::INPUT,},
  {"Max 1 (V)", "Neuron 1 maximum membrane potential (in V)", DefaultGUIModel::INPUT,},
  {"Min 2 (V)", "Neuron 2 minimum membrane potential (in V)", DefaultGUIModel::INPUT,},
  {"Max 2 (V)", "Neuron 2 maximum membrane potential (in V)", DefaultGUIModel::INPUT,},

  {"Scale 1-2", "Scale from neuron 1 to neuron 2", DefaultGUIModel::OUTPUT, },
  {"Offset 1-2", "Offset from neuron 1 to neuron 2", DefaultGUIModel::OUTPUT, },
  {"Scale 2-1", "Scale from neuron 2 to neuron 1", DefaultGUIModel::OUTPUT, },
  {"Offset 2-1", "Offset from neuron 2 to neuron 1", DefaultGUIModel::OUTPUT, },

  {"Scale 1-2", "Scale from neuron 1 to neuron 2", DefaultGUIModel::STATE,},
  {"Offset 1-2", "Offset from neuron 1 to neuron 2", DefaultGUIModel::STATE,},
  {"Scale 2-1", "Scale from neuron 2 to neuron 1", DefaultGUIModel::STATE,},
  {"Offset 2-1", "Offset from neuron 2 to neuron 1", DefaultGUIModel::STATE,},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

RTHybridAmplitudeScaleOffset::RTHybridAmplitudeScaleOffset(void)
  : DefaultGUIModel("RTHybridAmplitudeScaleOffset", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>RTHybrid Amplitude Scale Offset:</b><br>Given two neurons membrane potential minimum and maximum values, this module calculates the amplitude scale factor and offset between them, in both directions. i.e. if Neuron 1 membrane potential is multiplied by Scale 1-2 and added Offset 1-2, the result will be in the same amplitude range than Neuron 2 membrane potential, and viceversa.</p>");
  DefaultGUIModel::createGUI(vars,
                             num_vars); // this is required to create the GUI
  initParameters();
  update(INIT); // this is optional, you may place initialization code directly
                // into the constructor
  refresh();    // this is required to update the GUI with parameter and state
                // values
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

RTHybridAmplitudeScaleOffset::~RTHybridAmplitudeScaleOffset(void)
{
}

void
RTHybridAmplitudeScaleOffset::execute(void)
{
  double min1 = input(0);
  double max1 = input(1);
  double min2 = input(2);
  double max2 = input(3);

  double range1 = max1 - min1;
  double range2 = max2 - min2;

  s12 = range2 / range1;
  s21 = range1 / range2;

  o12 = min2 - (min1 * s12);
  o21 = min1 - (min2 * s21);

  output(0) = s12;
  output(1) = o12;
  output(2) = s21;
  output(3) = o21;

  return;
}

void
RTHybridAmplitudeScaleOffset::initParameters(void)
{
  s12 = 1.0;
  s21 = 1.0;
  o12 = 0.0;
  o21 = 0.0;
}

void
RTHybridAmplitudeScaleOffset::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag) {
    case INIT:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      setState("Scale 1-2", s12);
      setState("Offset 1-2", o12);
      setState("Scale 2-1", s21);
      setState("Offset 2-1", o21);
      break;

    case MODIFY:
      break;

    case UNPAUSE:
      break;

    case PAUSE:
      break;

    case PERIOD:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      break;

    default:
      break;
  }
}
