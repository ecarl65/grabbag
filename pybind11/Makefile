CXX=c++
CXXFLAGS=-O3 -Wall -shared -std=c++11 -fPIC  $(shell python3 -m pybind11 --includes)

example: example.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@$(shell python3-config --extension-suffix)

clean:
	rm -f example$(shell python3-config --extension-suffix)

