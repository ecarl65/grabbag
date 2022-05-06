/*
 * =====================================================================================
 *
 *       Filename:  example.cpp
 *
 *    Description:  Example pybind11
 *
 *        Version:  1.0
 *        Created:  04/22/2022 11:09:43 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i=1, int j=2) {
    return i + j;
}

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add1", &add, "A function that adds two numbers",
        py::arg("i") = 1, py::arg("j") = 2);

    using namespace pybind11::literals;
    m.def("add2", &add, "A function that adds two numbers", "i"_a = 1, "j"_a = 2);
}
