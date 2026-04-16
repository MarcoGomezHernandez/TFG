# RTXI plugin Makefile for the bidirectional chemical synapse BO module.
# Only build settings are defined here; RTXI provides the compilation rules.

PLUGIN_NAME = bidirectional_chemical_synapse_BO

# Public headers installed/used by the plugin.
HEADERS = bidirectional_chemical_synapse_BO.h\
          utils.hpp\

# Sources (includes BO + evaluation helpers).
SOURCES = bidirectional_chemical_synapse_BO.cpp\
          moc_bidirectional_chemical_synapse_BO.cpp\
          BO.cpp\
          evaluation.cpp\

# External dependencies (KFR DSP, Boost, NLOpt, oneTBB).
LIBS = -lkfr_dsp_avx2 -lboost_serialization -lboost_filesystem -lboost_system -lnlopt

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile

CXXFLAGS := $(filter-out -std=c++11,$(CXXFLAGS)) -std=c++20 -O2 -march=native -DNDEBUG -DUSE_NLOPT

LIMBO_DIR ?= ../Limbo/limbo
CXXFLAGS += -I$(LIMBO_DIR)/src -I/usr/include/eigen3 -I$(HOME)/.local/include
