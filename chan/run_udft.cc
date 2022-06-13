/*
 * =====================================================================================
 *
 *       Filename:  udft.cc
 *
 *    Description:  Testing overlap/save with FFTW in C to get the basic right before
 *                  using in the channelizer.
 *
 *        Version:  1.0
 *        Created:  05/23/2022 07:32:05 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eric Carlsen (carlsen.eric@gmail.com), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <iostream>
#include <complex>
#include <string>
#include <chrono>
#include <exception>
#include <getopt.h>
#include <fftw3.h>
#include "udft.hh"

using namespace std::chrono;

// {{{ main
int main(int argc, char **argv) {
  // Parse arguments, start with defaults.
  int downsamp = 8;
  int filt_ord = 8;
  int full_ord = 8;
  float samp_rate = 10e3;
  bool debug = false;
  bool write = false;
  int opt;
  while((opt = getopt(argc, argv, ":d:f::n::s::h::vw")) != -1) 
  { 
    switch(opt) 
    { 
      case 'h': 
        std::cout << "Run UDFT Channelizer.\n\n";
        std::cout << "Arguments:\n";
        std::cout << "\td - Downsample integer\n";
        std::cout << "\tf - Filter order (full length is order * downsample + 1)\n";
        std::cout << "\tn - Number of total samples (power of two, which is multiplied by downsamp)\n";
        std::cout << "\ts - Sample rate\n";
        std::cout << "\tv - Verbose output\n";
        std::cout << "\tw - Write output files\n";
        std::cout << "\th - This help\n";
        exit(EXIT_FAILURE);
        break;
      case 'd': 
        downsamp = atoi(optarg);
        break;
      case 'v': 
        debug = true;
        break;
      case 'w': 
        write = true;
        break;
      case 'f': 
        filt_ord = atoi(optarg);
        break; 
      case 'n': 
        full_ord = atoi(optarg);
        break; 
      case 's': 
        samp_rate = atof(optarg);
        break; 
      case ':': 
             std::cout << "option needs a value\n"; 
             break; 
      case '?': 
             std::cout << "unknown option: " << optopt << std::endl;
             break; 
    } 
  } 
  int n_full = pow(2, full_ord) * downsamp;  // Total number of samples
  int n_filt = filt_ord * downsamp + 1;

  // Set up threads
  int init_threads = fftw_init_threads();
  if (init_threads == 0) throw std::runtime_error("Could not init FFTW threads\n");
  fftw_plan_with_nthreads(downsamp);

  // Set up channelizer
  UDFT channelizer(downsamp, n_filt, samp_rate, write, debug);

  // Set up input signal
  float *full_in = fftwf_alloc_real(n_full);
  float chirp_period = n_full / samp_rate / 2; // This many *full* cycles of chirp
  make_chirp(full_in, n_full, samp_rate, chirp_period);

  // Run and time result
  auto start = high_resolution_clock::now();
  auto full_out = channelizer.run(full_in, n_full);
  auto duration = duration_cast<nanoseconds>(high_resolution_clock::now() - start);
  std::cout << "Elapsed time for run call: " << duration.count() * 1e-9 << " seconds" << std::endl;

  // Free input memory
  fftwf_free(full_in);
}
// }}}
