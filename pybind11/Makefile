CXX=c++
CXXFLAGS=-O3 -Wall -fopenmp -std=c++11 -fPIC  $(shell pybind11-config --includes)
EXT=$(shell python3.8-config --extension-suffix)

.PHONY: clean

all: example iteration_mod pet

example: example.cpp
	$(CXX) -shared $(CXXFLAGS) $(CPPFLAGS) $< -o $@$(EXT)

iteration_mod: iteration_mod.cpp
	$(CXX) -shared $(CXXFLAGS) $< -o $@$(EXT)

pet: pet.o Pet.o
	$(CXX) -shared -o $@$(EXT) $(CXXFLAGS) $^

pet.o: pet.cpp
	$(CXX) $(CXXFLAGS) -c $<

Pet.o: Pet.cpp Pet.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.so *.o example

