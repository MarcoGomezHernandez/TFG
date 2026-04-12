PLUGIN_NAME = bidirectional_chemical_synapse_genetic

HEADERS = bidirectional_chemical_synapse_genetic.h\
          utils.hpp\

SOURCES = bidirectional_chemical_synapse_genetic.cpp\
          moc_bidirectional_chemical_synapse_genetic.cpp\
          bayesian.cpp\
          fitness.cpp\

LIBS = -lkfr_dsp_avx2 -lboost_serialization -lboost_filesystem -lboost_system -lnlopt -ltbb

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile

CXXFLAGS := $(filter-out -std=c++11,$(CXXFLAGS)) -std=c++20 -O2 -march=native -DNDEBUG -DUSE_NLOPT -DUSE_TBB -DUSE_TBB_ONEAPI

LIMBO_DIR ?= ../Limbo/limbo
CXXFLAGS += -I$(LIMBO_DIR)/src -I/usr/include/eigen3 -I$(HOME)/.local/include
