CXX = g++
CXXFLAGS = -Wall -O2 -I/usr/local/Neun/0.4.0 -std=c++20
LDLIBS = -lyaml-cpp

TARGET = programa
SRC = parametrize_chemicalSynapsis_PD-HR.cpp scaling.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDLIBS)

clean:
	rm -f $(TARGET) *.o
