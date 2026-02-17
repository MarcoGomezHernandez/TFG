CXX = g++
CXXFLAGS = -Wall -O2 -I/usr/local/Neun/0.4.0 -std=c++20

# Targets
CONSTS_CALCULATOR = consts_calculator
GENETIC = genetic

# Source files
CONSTS_SRCS = consts_calculator.cpp scaling.cpp
GENETIC_SRCS = main.cpp genetic.cpp fitness.cpp scaling.cpp

# Object files
CONSTS_OBJS = $(CONSTS_SRCS:.cpp=.o)
GENETIC_OBJS = $(GENETIC_SRCS:.cpp=.o)

# Default target: build both executables
all: $(CONSTS_CALCULATOR) $(GENETIC)

# Build consts_calculator
$(CONSTS_CALCULATOR): $(CONSTS_OBJS)
	$(CXX) $(CXXFLAGS) $(CONSTS_OBJS) -o $(CONSTS_CALCULATOR) $(LDLIBS)

# Build genetic
$(GENETIC): $(GENETIC_OBJS)
	$(CXX) $(CXXFLAGS) $(GENETIC_OBJS) -o $(GENETIC) $(LDLIBS)

# Pattern rule for object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean specific target
clean-consts:
	rm -f $(CONSTS_CALCULATOR) $(CONSTS_OBJS)

clean-genetic:
	rm -f $(GENETIC) $(GENETIC_OBJS)

# Clean everything
clean:
	rm -f $(CONSTS_CALCULATOR) $(GENETIC) *.o

# Phony targets
.PHONY: all clean clean-consts clean-genetic
