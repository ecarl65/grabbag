/*
 * =====================================================================================
 *
 *       Filename:  Pet.h
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
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Python.h>
namespace py = pybind11;

struct Pet {
  Pet(const std::string &name);
  void setName(const std::string &name_);
  const std::string &getName();
  py::array nump(py::array arr);


  std::string name;
};

