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
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <tgmath.h>
#include <fftw3.h>
#include <string.h>
#include "udft.hh"

using namespace std::chrono;

// {{{ main
int main(int argc, char **argv) {
  // Values on which all others depend
  int downsamp = 8; // Downsample factor
  int filt_ord = 8;
  int full_ord = 256;
  float samp_rate = 10e3;
  bool debug = false;
  bool write = false;

  int opt;

  // put ':' in the starting of the string so that program can distinguish between '?' and ':' 
  while((opt = getopt(argc, argv, ":d:f::n::s::h::vw")) != -1) 
  { 
    switch(opt) 
    { 
      case 'h': 
        printf("Run UDFT Channelizer.\n\n");
        printf("Arguments:\n");
        printf("\td - Downsample integer\n");
        printf("\tf - Filter order (full length is order * downsample + 1)\n");
        printf("\tn - Number of total samples (multiple of downsamp)\n");
        printf("\ts - Sample rate\n");
        printf("\tv - Verbose output\n");
        printf("\tw - Write output files\n");
        printf("\th - This help\n");
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
             printf("option needs a value\n"); 
             break; 
      case '?': 
             printf("unknown option: %c\n", optopt);
             break; 
    } 
  } 
  int n_full = full_ord * downsamp;  // Total number of samples
  int n_filt = filt_ord * downsamp + 1;

  // Set up threads
  int init_threads = fftw_init_threads();
  if (init_threads == 0) exit(EXIT_FAILURE);
  fftw_plan_with_nthreads(downsamp);

  UDFT channelizer(downsamp, n_full, n_filt, samp_rate, write, debug);
  auto start = high_resolution_clock::now();
  channelizer.run();
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  std::cout << "Elapsed time for run call: " << duration.count() * 1e-6 << " seconds" << std::endl;

}
// }}}
