
CXX=g++
CXXFLAGS=-g -Ofast -std=c++11 -Wall -fopenmp -fPIC $(shell pybind11-config --includes)
LDFLAGS=-lm -lfftw3f -lfftw3 -lfftw3_omp -lfftw3f_omp $(shell python3.8-config --ldflags)
EXT=$(shell python3.8-config --extension-suffix)

all: run_udft pyudft

pyudft: pyudft.cc udft.o
	$(CXX) -shared -o $@$(EXT) $(CXXFLAGS) $^ $(LDFLAGS)

run_udft: run_udft.o udft.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run_udft.o: run_udft.cc udft.o
	$(CXX) $(CXXFLAGS) -c $<

udft.o: udft.cc udft.hh utils.hh
	$(CXX) $(CXXFLAGS) -c $<

.phony: clean

clean:
	rm -f *.bin *.o run_udft

