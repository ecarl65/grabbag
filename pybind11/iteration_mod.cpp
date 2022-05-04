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
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Python.h>
namespace py = pybind11;
// py::module nn = py::module::import("iteration");


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

    m.def("nump", &nump,
      "the function which loops over a numpy array");
}
