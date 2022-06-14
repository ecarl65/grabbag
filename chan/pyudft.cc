/*
 * =====================================================================================
 *
 *       Filename:  pyudft.cpp
 *
 *    Description:  Object oriented version of pybind11
 *
 *        Version:  1.0
 *        Created:  05/06/2022 08:49:35 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <string>
#include <pybind11/pybind11.h>
#include "pudft.hh"

namespace py = pybind11;

PYBIND11_MODULE(pyudft, m) {
    py::class_<PUDFT>(m, "pudft")
        .def(py::init<int, int, float, bool, bool>(), "Args: Downsample, Filter Size, Sample Rate, Write Binary Output, Verbose")
        .def("poly_filt_design", &PUDFT::poly_filt_design, "Design and transform the filter")
        .def("run", &PUDFT::run, "Run with python array input")
        .def("pyrun", &PUDFT::pyrun, "Run with C array input")
        .def("get_filter", &PUDFT::get_filter, "Get the designed filter");
}

