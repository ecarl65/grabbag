CXX=g++
NVCC=nvcc
TARGETS=add add_cuda add_cuda1 add_cuda2

all: $(TARGETS)

add: add.cpp
	$(CXX) -o $@ $^

add_cuda: add_cuda.cu
	$(NVCC) -o $@ $^

add_cuda1: add_cuda1.cu
	$(NVCC) -o $@ $^

add_cuda2: add_cuda2.cu
	$(NVCC) -o $@ $^
clean:
	rm -f $(TARGETS)
