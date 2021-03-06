/*
 * =====================================================================================
 *
 *       Filename:  udft.hh
 *
 *    Description:  Header file for Floating Point UDFT
 *
 *        Version:  1.0
 *        Created:  06/12/2022 03:13:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <iostream>
#include <complex>
#include <vector>
#include <fftw3.h>
#include <stdexcept>
#include "pudft.hh"

// {{{ pyrun(py::array)
py::array PUDFT::pyrun(py::array indata)
{

  auto indata_obj_prop = indata.request();

  //initialize values
  float *in_vals = (float*) indata_obj_prop.ptr;

  int n_full = static_cast<int>(indata_obj_prop.shape[0]);

  auto output = run(in_vals, n_full);

  py::array ret =  py::cast(output); //py::array(vect_arr.size(), vect_arr.data());
  return ret;

}
// }}}


// {{{ get_filter
py::array PUDFT::get_filter()
{

  std::vector<float> filter(n_filt);
  for (int m = 0; m < n_filt; m++) filter[m] = filt[m];

  py::array ret =  py::cast(filter); //py::array(vect_arr.size(), vect_arr.data());
  return ret;

}
// }}}
