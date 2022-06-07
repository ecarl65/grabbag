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

  // FFT parameters
  struct fft_config fwd_c = {
    {ncols}, 1, M, 1, ncols_fft, M, 1, NULL, NULL
  };
  struct fft_config filt_c = {
    {ncols}, 1, M, ncols, ncols_fft, 1, 1, NULL, NULL
  };
  struct fft_config inv_c = {
    {ncols}, 1, M, ncols_fft, ncols, 1, 1, NULL, NULL
  };
  struct fft_config col_c = {
    {M}, 1, ncols, 1, 1, ncols, ncols, NULL, NULL
  };

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

  // Make input chirp
  for (int m = 0; m < Nfull; m++) {
    /* full_in[m] = cos(2 * pi * fc * m * Ts) + coarse_gaussian() / 10.0; */
    full_in[m] = sin(2 * pi * (2 * Fs / Tfull * pow(Ts * m, 2) / 2)) + coarse_gaussian() / 10.0;
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
  fftw_plan psig = fftw_plan_many_dft_r2c(fwd_c.rank, fwd_c.n_size, fwd_c.howmany, full_in,
                                          fwd_c.inembed, fwd_c.istride, fwd_c.idist, fft_in,
                                          fwd_c.onembed, fwd_c.ostride, fwd_c.odist, FFTW_ESTIMATE);

  // Filter FFT
  fftw_plan pfilt = fftw_plan_many_dft_r2c(filt_c.rank, filt_c.n_size, filt_c.howmany, filt_full,
                                           filt_c.inembed, filt_c.istride, filt_c.idist, fft_filt,
                                           filt_c.onembed, filt_c.ostride, filt_c.odist, FFTW_ESTIMATE);
  fftw_execute(pfilt);
  
  // IFFT plan
  fftw_plan pinv = fftw_plan_many_dft_c2r(inv_c.rank, inv_c.n_size, inv_c.howmany, fft_mult,
                                          inv_c.inembed, inv_c.istride, inv_c.idist, conv_out,
                                          inv_c.onembed, inv_c.ostride, inv_c.odist, FFTW_ESTIMATE);

  // TODO: this section will be in the loop

  // Forward FFT
  fftw_execute(psig);

  // Compute circular convolution
  for (int m = 0; m < Nfft_h; m++) {
    fft_mult[m] = fft_filt[m] * fft_in[m];
  }
  fftw_execute(pinv);

  // Perform DFT down columns to get channelized output
  fftw_plan pudft = fftw_plan_many_dft_r2c(col_c.rank, col_c.n_size, col_c.howmany, conv_out,
                                           col_c.inembed, col_c.istride, col_c.idist, udft,
                                           col_c.onembed, col_c.ostride, col_c.odist, FFTW_ESTIMATE);
  fftw_execute(pudft);

  // End loop

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



