CXX = g++
CXXFLAGS = -O3 -std=c++20
PYTHON = python

BIN_DIR = bin
SRC_DIR = src
SCRIPTS_DIR = scripts
DATA_DIR = data
FIGURES_DIR = figures

HRM_TARGET = hindmarsh_rose_improved
HRM_SIM_BIN = $(BIN_DIR)/$(HRM_TARGET)_sim
HRM_SIM_SRC = $(SRC_DIR)/$(HRM_TARGET)_sim.cpp
HRM_PLOT_SCRIPT = $(SCRIPTS_DIR)/$(HRM_TARGET)_plot.py
HRM_DATA_FILE = $(DATA_DIR)/$(HRM_TARGET)_data.csv
HRM_FIGURE_FILE = $(FIGURES_DIR)/$(HRM_TARGET)_fig.png

TARGETS = $(HRM_TARGET)

all: $(TARGETS)

$(HRM_TARGET): $(HRM_SIM_BIN) $(HRM_PLOT_SCRIPT) | $(DATA_DIR) $(FIGURES_DIR)
	./$(HRM_SIM_BIN)
	$(PYTHON) $(HRM_PLOT_SCRIPT)

$(HRM_SIM_BIN): CXXFLAGS += -I/usr/local/Neun/0.4.0
$(HRM_SIM_BIN): $(HRM_SIM_SRC) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BIN_DIR) $(DATA_DIR) $(FIGURES_DIR):
	mkdir -p $@

clean:
	rm -f $(HRM_SIM_BIN) $(HRM_DATA_FILE) $(HRM_FIGURE_FILE)

.PHONY: all $(TARGETS) clean
