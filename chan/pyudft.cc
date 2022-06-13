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
#include "udft.hh"

namespace py = pybind11;

PYBIND11_MODULE(pyudft, m) {
    py::class_<UDFT>(m, "udft")
        .def(py::init<int, int, float, bool, bool>())
        .def("poly_filt_design", &UDFT::poly_filt_design, "Design and transform the filter")
        .def("run", static_cast<void (UDFT::*)(py::array)>(&UDFT::run), "Run with python array input")
        .def("run", static_cast<void (UDFT::*)(float *, int)>(&UDFT::run), "Run with C array input");
}

