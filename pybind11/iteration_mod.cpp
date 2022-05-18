/*
 * =====================================================================================
 *
 *       Filename:  vec.cc
 *
 *    Description:  Example vector/numpy array form stackoverflow (29956898)
 *
 *        Version:  1.0
 *        Created:  05/04/2022 07:24:15 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Python.h>
namespace py = pybind11;

py::array_t<double> add_arrays(py::array_t<double> input1, py::array_t<double> input2) {
    py::buffer_info buf1 = input1.request(), buf2 = input2.request();

    if (buf1.ndim != 1 || buf2.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    if (buf1.size != buf2.size)
        throw std::runtime_error("Input shapes must match");

    /*  No pointer is passed, so NumPy will allocate the buffer */
    auto result = py::array_t<double>(buf1.size);

    py::buffer_info buf3 = result.request();

    double *ptr1 = static_cast<double *>(buf1.ptr);
    double *ptr2 = static_cast<double *>(buf2.ptr);
    double *ptr3 = static_cast<double *>(buf3.ptr);

    for (ssize_t idx = 0; idx < buf1.shape[0]; idx++)
        ptr3[idx] = ptr1[idx] + ptr2[idx];

    return result;
}

py::array nump(py::array arr){

    auto arr_obj_prop = arr.request();
    //initialize values
    double *vals = (double*) arr_obj_prop.ptr;

    unsigned int shape_1 = arr_obj_prop.shape[0];
    unsigned int shape_2 = arr_obj_prop.shape[1];


    std::vector<std::vector <double>> vect_arr( shape_1, std::vector<double> (shape_2));

    for(unsigned int i = 0; i < shape_1; i++){
      for(unsigned int j = 0; j < shape_2; j++){
        vect_arr[i][j] = vals[i*shape_2 + j] * 2;
      }
    }   

    py::array ret =  py::cast(vect_arr); //py::array(vect_arr.size(), vect_arr.data());
    return ret;

}

PYBIND11_MODULE(iteration_mod, m) {

    m.doc() = "pybind11 module for iterating over generations";

    m.def("nump", &nump, "the function which loops over a numpy array");
    m.def("add_arrays", &add_arrays, "Add two NumPy arrays");
}
