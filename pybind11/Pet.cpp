/*
 * =====================================================================================
 *
 *       Filename:  Pet.cpp
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
#include <vector>
#include "Pet.h"

Pet::Pet(const std::string &name) : name(name) { }
void Pet::setName(const std::string &name_) { 
  name = name_;
}
const std::string &Pet::getName() {
  return name;
}
py::array Pet::nump(py::array arr) {

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
