CXX=c++
# CXXFLAGS=-O3 -Wall -fopenmp -shared -std=c++11 -fPIC  $(shell python3 -m pybind11 --includes)
CXXFLAGS=-O3 -Wall -fopenmp -shared -std=c++11 -fPIC  $(shell pybind11-config --includes)
SUFFIX=$(shell python3.8-config --extension-suffix)

.PHONY: clean

all: example iteration_mod pet

example: example.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@$(shell python3.8-config --extension-suffix)

iteration_mod: iteration_mod.cpp
	$(CXX) $(CXXFLAGS) $< -o $@$(shell python3.8-config --extension-suffix)

pet: Pet.cpp
	$(CXX) $(CXXFLAGS) $< -o $@$(shell python3.8-config --extension-suffix)

clean:
	rm -f example$(SUFFIX) iteration_mod$(SUFFIX) pet$(SUFFIX)

