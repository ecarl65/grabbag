/*
 * =====================================================================================
 *
 *       Filename:  Pet_obj.h
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

struct Pet {
  Pet(const std::string &name);
  void setName(const std::string &name_);
  const std::string &getName();

  std::string name;
};

