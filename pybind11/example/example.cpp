/*
 * =====================================================================================
 *
 *       Filename:  example.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/21/2022 07:45:04 PM
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

int add(int i, int j) {
      return i + j;
}

PYBIND11_MODULE(example, m) {
      m.doc() = "pybind11 example plugin"; // optional module docstring

      m.def("add1", &add, "A function that adds two numbers",
          py::arg("i")  , py::arg("j"));
      using namespace pybind11::literals;
      m.def("add2", &add, "A function that adds two numbers", "i"_a, "j"_a);
          
}

