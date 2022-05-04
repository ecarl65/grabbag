CXX=c++
CXXFLAGS=-O3 -Wall -fopenmp -shared -std=c++11 -fPIC  $(shell python3 -m pybind11 --includes)

.PHONY: clean

all: example iteration_mod

example: example.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@$(shell python3-config --extension-suffix)

iteration_mod: iteration_mod.cpp
	$(CXX) $(CXXFLAGS) $< -o $@$(shell python3-config --extension-suffix)

clean:
	rm -f example$(shell python3-config --extension-suffix) iteration_mod$(shell python3-config --extension-suffix)

