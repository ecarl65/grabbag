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
#include "Pet.h"

Pet::Pet(const std::string &name) : name(name) { }
void Pet::setName(const std::string &name_) { 
  name = name_;
}
const std::string &Pet::getName() {
  return name;
}

