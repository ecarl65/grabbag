/*
 * =====================================================================================
 *
 *       Filename:  osconv.c
 *
 *    Description:  Testing overlap/save with FFTW in C to get the basic right before
 *                  using in the channelizer.
 *
 *        Version:  1.0
 *        Created:  05/23/2022 07:32:05 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>
#include <tgmath.h>
#include <string.h>
#include "utils.h"

int main(int argc, char **argv) {
  // Variables
  const int M = 8; // Downsample factor
  const int Nfull = 256 * M;  // Total number of samples
  const int Nds = Nfull / M;
  const int Nfilt = 16 * M + 1;
  const int Nfilt_half = (Nfilt - 1) >> 1;
  const int Noutdelay = Nfilt_half / M;
  const int Nfill = Nfull - Nfilt + 1;
  const double pi = acosf(-1);
  const double Fs = 10e3;
  const double Ts = 1.0 / Fs;
  const double Tfull = Nfull * Ts;
  const double fc = 2600;
  const double fcutoff = Fs / (2 * M);
  const int nrows = M;
  const int ncols = Nds;
  const int ncols_fft = ((Nfull / M) / 2 + 1);
  const int nrows_fft = M / 2 + 1;
  const int Nfft_h = ncols_fft * M;
  const int Nfft_v = nrows_fft * ncols;
  srand(time(NULL));

  // Arrays
  double *full_in = (double *) fftw_malloc(sizeof(double) * Nfull);
  double filt[Nfilt];
  double *filt_full = (double *) fftw_malloc(sizeof(double) * Nfull);
  fftw_complex *fft_in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nfft_h);
  fftw_complex *fft_filt = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nfft_h);
  double *conv_out = (double *) fftw_malloc(sizeof(double) * Nfull);
  fftw_complex *fft_mult = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nfft_h);
  fftw_complex *udft = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nfft_v);

  // Make filter
  int outidx = 0;
  double filt_sum = 0;
  for (int m = -Nfilt_half; m <= Nfilt_half; m++) {
    // Sinc low pass filter at twice the SOI freq
    if (m == 0) {
      filt[outidx] = 2 * fcutoff / Fs;
    } else {
      filt[outidx] = sin(2 * pi * fcutoff * m / Fs) / (m * pi);
    }

    // Hamming window and normalization
    filt[outidx] *= (0.54 - 0.46 * cos(2 * pi * outidx / Nfilt));
    filt_sum += filt[outidx++];
  }
  for (int m = 0; m < Nfilt; m++) filt[m] /= filt_sum;  // Normalize the filter

  // Make input sinusoid
  for (int m = 0; m < Nfull; m++) {
    if (m < Nfill) {
      /* full_in[m] = cos(2 * pi * fc * m * Ts) + coarse_gaussian() / 10.0; */
      full_in[m] = sin(2 * pi * (2 * Fs / Tfull * pow(Ts * m, 2) / 2)) + coarse_gaussian() / 10.0;
    } else {
      full_in[m] = 0;
    }
  }

  // Filter polyphase decomposition
  for (int n = 0; n < ncols; n++) {
    for (int rho = 0; rho < M; rho++) {
      int inidx = n * M - rho;
      int outidx = rho * ncols + n;
      if (inidx < 0 || inidx >= Nfilt) {
        filt_full[outidx] = 0;
      } else {
        filt_full[outidx] = filt[inidx];
      }
      outidx++;
    }
  }

  // Data polyphase FFT
  int n_size[] = {ncols};
  int rank = 1;
  int howmany = M;
  int idist = 1;
  int odist = ncols_fft;
  int istride = M;
  int ostride = 1;
  int *inembed = NULL;
  int *onembed = NULL;
  fftw_plan psig = fftw_plan_many_dft_r2c(rank, n_size, howmany, full_in, inembed, istride, idist, fft_in,
                                            onembed, ostride, odist, FFTW_ESTIMATE);

  // Filter FFT
  n_size[0] = ncols;
  rank = 1;
  howmany = M;
  idist = ncols;
  odist = ncols_fft;
  istride = 1;
  ostride = 1;
  inembed = NULL;
  onembed = NULL;
  fftw_plan pfilt = fftw_plan_many_dft_r2c(rank, n_size, howmany, filt_full, inembed, istride, idist, fft_filt,
                                             onembed, ostride, odist, FFTW_ESTIMATE);
  fftw_execute(pfilt);
  
  // IFFT plan
  n_size[0] = ncols;  // LOGICAL size, not PHYSICAL size of output array
  rank = 1;
  howmany = M;
  idist = ncols_fft;
  odist = ncols;
  istride = 1;
  ostride = 1;
  inembed = NULL;
  onembed = NULL;
  fftw_plan pinv = fftw_plan_many_dft_c2r(rank, n_size, howmany, fft_mult, inembed, istride, idist, conv_out,
                                            onembed, ostride, odist, FFTW_ESTIMATE);

  // Forward FFT
  fftw_execute(psig);

  // Compute circular convolution
  for (int m = 0; m < Nfft_h; m++) {
    fft_mult[m] = fft_filt[m] * fft_in[m];
  }
  fftw_execute(pinv);

  // Perform DFT down columns to get channelized output
  n_size[0] = M;
  rank = 1;
  howmany = ncols;
  idist = 1;
  odist = 1;
  istride = ncols;
  ostride = ncols;
  inembed = NULL;
  onembed = NULL;
  fftw_plan pudft = fftw_plan_many_dft_r2c(rank, n_size, howmany, conv_out, inembed, istride, idist, udft,
                                             onembed, ostride, odist, FFTW_ESTIMATE);
  fftw_execute(pudft);

  // Write outputs
  write_out("filtered.bin", (void *) &conv_out[0], sizeof(double), Nfull);
  write_out("input.bin", (void *) &full_in[0], sizeof(double), Nfull);
  write_out("filter.bin", (void *) &filt[0], sizeof(double), Nfilt);
  write_out("channelized.bin", (void *) &udft[0], sizeof(fftw_complex), Nfft_v);
  write_out("fftdata.bin", (void *) &fft_in[0], sizeof(fftw_complex), Nfft_h);
  write_out("fftfilt.bin", (void *) &fft_filt[0], sizeof(fftw_complex), Nfft_h);

  // Free memory
  fftw_free(full_in);
  fftw_free(filt_full);
  fftw_free(udft);
  fftw_free(fft_in);
  fftw_free(fft_filt);
  fftw_free(conv_out);
  fftw_free(fft_mult);
  fftw_destroy_plan(pfilt);
  fftw_destroy_plan(psig);
  fftw_destroy_plan(pinv);
  fftw_destroy_plan(pudft);
  fftw_cleanup();
}



