CXX=c++
# CXXFLAGS=-O3 -Wall -shared -std=c++11 -fPIC -I/usr/include/python3.6m -I/home/ecarlsen/pvenv/lib64/python3.6/site-packages/pybind11/include
CXXFLAGS=-O3 -Wall -shared -std=c++11 -fPIC  $(shell python3 -m pybind11 --includes)

example: example.cpp
	# $(CXX) -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) example.cpp -o example$(python3-config --extension-suffix)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@$(shell python3-config --extension-suffix)

