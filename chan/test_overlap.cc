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
#include <fftw3.h>
#include <stdexcept>
#include "utils.hh"

using cfloat = std::complex<float>;

int main(int argc, char **argv) {
  int M = 4;
  int K = 8;
  int I = 2;
  int N = 256;
  int n_fft = N / M;
  int n_out_fft = n_fft / 2 + 1;
  int extra = (I - 1) * M;
  int n_tot = n_out_fft * K;

  float *buffer_in = fftwf_alloc_real(n_fft + extra);
  cfloat *fft_in = reinterpret_cast<cfloat*>(fftwf_alloc_complex(n_tot));

  // {{{ FFT Plan Setup
  // Different approach, not upsampling/replicating
  int n_size[] = { n_fft };
  int rank = 1;
  int howmany = K;
  int idist = 1;
  int odist = n_out_fft;
  int istride = M;
  int ostride = 1;
  int *inembed = NULL;
  int *onembed = NULL;

  // Data polyphase FFT
  fftwf_plan psig = fftwf_plan_many_dft_r2c(rank, n_size, howmany, buffer_in,
      inembed, istride, idist,
      reinterpret_cast<fftwf_complex*>(fft_in),
      onembed, ostride, odist, FFTW_MEASURE);

  // Fill in input data
  for (size_t m = 0; m < n_fft + extra; m++) {
    buffer_in[m] = sin(2 * 3.14159 * 0.3 * m);
  }

  // Do FFT
  fftwf_execute(psig);

  // Cleanup
  fftwf_free(buffer_in);
  fftwf_free(fft_in);
  fftwf_destroy_plan(psig);
}
// }}}



