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

#include "rthybrid_electrical_synapse.h"
#include <iostream>
#include <main_window.h>

#define SM_ELECTRICAL_N_G 1
#define SM_ELECTRICAL_G 0

extern "C" Plugin::Object*
createRTXIPlugin(void)
{
  return new RTHybridElectricalSynapse();
}

static DefaultGUIModel::variable_t vars[] = {
    { "Current (nA)", "Synaptic current (in nA)", DefaultGUIModel::OUTPUT, },

    {"g (uS)", "Conductance (uS)", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

    { "Post-synaptic Voltage (V)", "", DefaultGUIModel::INPUT,},
    { "Pre-synaptic Voltage (V)", "", DefaultGUIModel::INPUT,},

    { "Scale (Pre to Post)", "Scale from the pre-synaptic neuron to the post-synaptic one, as given by the RTHybridAmplitudeScale module.", DefaultGUIModel::INPUT,},
    { "Offset (Pre to Post)", "Offset from the pre-synaptic neuron to the post-synaptic one, as given by the RTHybridAmplitudeScale module.", DefaultGUIModel::INPUT,},
   
    { "Current (nA)", "", DefaultGUIModel::STATE,},
    { "offset", "", DefaultGUIModel::STATE,},
    { "scale", "", DefaultGUIModel::STATE,},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

RTHybridElectricalSynapse::RTHybridElectricalSynapse(void)
  : DefaultGUIModel("RTHybrid electrical synapse model", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>RTHybrid electrical synapse model</b><br>RTHybrid module for RTXI that implements a gap junction electrical synapse model.</p>");
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

RTHybridElectricalSynapse::~RTHybridElectricalSynapse(void)
{
}

void RTHybridElectricalSynapse::sm_electrical (double v_post, double v_pre, double * ret) {
    v_pre = v_pre * scale + offset;

    *ret = (g[SM_ELECTRICAL_G] * (v_post - v_pre));

    return;
}

void
RTHybridElectricalSynapse::execute(void)
{
    scale = input(2);
    offset = input(3) * 1000;

    if (scale == 0) {
        scale = 1;
        offset = 0;
    }

    sm_electrical(input(0) * 1000, input(1) * 1000, &i);

    output(0) = i;

    return;
}

void
RTHybridElectricalSynapse::initParameters(void)
{
    g[SM_ELECTRICAL_G] = 0;

    scale = 0;
    offset = 0;
}

void
RTHybridElectricalSynapse::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag) {
    case INIT:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      setParameter("g (uS)", g[SM_ELECTRICAL_G]);
      
      setState("Current (nA)", i);
      setState("offset", offset);
      setState("scale", scale);
      break;

    case MODIFY:
      g[SM_ELECTRICAL_G] = getParameter("g (uS)").toDouble();

      scale = 0;
      offset = 0;
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

void
RTHybridElectricalSynapse::customizeGUI(void)
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
RTHybridElectricalSynapse::aBttn_event(void)
{
}

void
RTHybridElectricalSynapse::bBttn_event(void)
{
}
