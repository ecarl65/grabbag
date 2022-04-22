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

#include <pybind11/pybind11.h>

int add(int i, int j) {
    return i + j;
}

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
}
