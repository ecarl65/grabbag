
CXX=g++
CXXFLAGS=-g -Ofast -std=c++11 -Wall -fopenmp -fPIC
PCXXFLAGS=$(shell pybind11-config --includes)
LDFLAGS=-lfftw3f -lfftw3 -lfftw3_omp -lfftw3f_omp -lm
EXT=$(shell python3-config --extension-suffix)

all: run_udft pyudft

pyudft: pyudft.o udft.o pudft.o
	$(CXX) -shared -o $@$(EXT) $(CXXFLAGS) $(PCXXFLAGS) $^ $(LDFLAGS) $(shell python3-config --ldflags)

pyudft.o: pyudft.cc
	$(CXX) $(CXXFLAGS) $(PCXXFLAGS) -c $<

run_udft: run_udft.o udft.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run_udft.o: run_udft.cc udft.o
	$(CXX) $(CXXFLAGS) -c $<

udft.o: udft.cc udft.hh utils.hh
	$(CXX) $(CXXFLAGS) -c $<

pudft.o: pudft.cc pudft.hh utils.hh
	$(CXX) $(CXXFLAGS) $(PCXXFLAGS) -c $<

.phony: clean

clean:
	rm -f *.bin *.o run_udft *.so

