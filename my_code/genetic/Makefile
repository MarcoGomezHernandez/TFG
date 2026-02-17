CXX = g++
CXXFLAGS = -Wall -O2 -I/usr/local/Neun/0.4.0 -std=c++20

# Targets
CONSTS_CALCULATOR = consts_calculator
GENETIC = genetic

# Source files
CONSTS_SRCS = consts_calculator.cpp scaling.cpp
GENETIC_SRCS = main.cpp scaling.cpp

# Default target: build both executables
all: $(CONSTS_CALCULATOR) $(GENETIC)

# Build consts_calculator directly from sources (no .o)
$(CONSTS_CALCULATOR): $(CONSTS_SRCS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

# Build genetic directly from sources (no .o)
$(GENETIC): $(GENETIC_SRCS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

# Clean specific target
clean-consts:
	rm -f $(CONSTS_CALCULATOR)

clean-genetic:
	rm -f $(GENETIC)

# Clean everything
clean:
	rm -f $(CONSTS_CALCULATOR) $(GENETIC)

# Phony targets
.PHONY: all clean clean-consts clean-genetic
