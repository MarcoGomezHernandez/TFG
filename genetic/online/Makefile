PLUGIN_NAME = bidirectional_chemical_synapse_genetic

HEADERS = bidirectional_chemical_synapse_genetic.h\
          fitness.hpp\

SOURCES = bidirectional_chemical_synapse_genetic.cpp\
          moc_bidirectional_chemical_synapse_genetic.cpp\
          genetic.cpp\
          fitness.cpp\

LIBS = -lkfr_dsp

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile

# Se ejecuta después de que RTXI inyecte sus banderas
CXXFLAGS := $(filter-out -std=c++11,$(CXXFLAGS)) -std=c++20
