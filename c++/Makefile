
CXX=g++
# CFLAGS=-g -DDEBUG
CFLAGS=-g -std=c++11
LIBS=-lm -lfftw3f -lfftw3 -lfftw3_omp -lfftw3f_omp

all: fudft

fudft: fudft.cc fudft.hh utils.hh
	$(CXX) -o $@ $< $(CFLAGS) $(LIBS)

.phony: clean

clean:
	rm -f *.bin fudft

