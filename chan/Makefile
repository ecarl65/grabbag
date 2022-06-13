
CXX=g++
CXXFLAGS=-g -std=c++11 -Wall
LDFLAGS=-lm -lfftw3f -lfftw3 -lfftw3_omp -lfftw3f_omp

all: run_udft

run_udft: run_udft.o udft.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

run_udft.o: run_udft.cc udft.o
	$(CXX) $(CXXFLAGS) -c $<

udft.o: udft.cc udft.hh utils.hh
	$(CXX) $(CXXFLAGS) -c $<

.phony: clean

clean:
	rm -f *.bin *.o run_udft

