CXX=c++
CXXFLAGS=-O3 -Wall -fopenmp -std=c++11 -fPIC  $(shell pybind11-config --includes)
EXT=$(shell python3.8-config --extension-suffix)

.PHONY: clean

all: example iteration_mod pet

example: example.cpp
	$(CXX) -shared $(CXXFLAGS) $(CPPFLAGS) $< -o $@$(shell python3.8-config --extension-suffix)

iteration_mod: iteration_mod.cpp
	$(CXX) -shared $(CXXFLAGS) $< -o $@$(shell python3.8-config --extension-suffix)

pet: pet.o Pet.o
	$(CXX) -shared -o $@$(shell python3.8-config --extension-suffix) $(CXXFLAGS) $< 

pet.o: pet.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

Pet.o: Pet.cpp Pet.h
	$(CXX) -O3 -Wall -fopenmp -std=c++11 -fPIC -c $< -o $@

clean:
	rm -f *.so *.o

