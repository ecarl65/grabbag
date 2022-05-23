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
#include <complex.h>
#include <fftw3.h>
#include <tgmath.h>

int main(int argc, char **argv) {
  // Make input sinusoid
  const int Nbuf = 1024;
  const int Nbufs = 20;
  const int Nfull = Nbuf * Nbufs;
  const double pi = acosf(-1);
  const double Fs = 10e3;
  const double Ts = 1.0 / Fs;
  double *full_in = (double *) fftw_malloc(sizeof(double) * Nfull);
  for (int m = 0; m < Nfull; m++) {
    full_in[m] = cos(2 * pi * 100 * m * Ts);
  }
  double *full_out = fftw_malloc(sizeof(double) * Nfull);

  // Make filter
  const int Nfilt = 128;
  double *filt = (double *) fftw_malloc(sizeof(double) * Nbuf);
  for (int m = 0; m < Nfilt; m++) {
    double sinc_idx = (m - 63.5) * 0.15 * 2 * pi;
    filt[m] = sin(sinc_idx) / sinc_idx;
  }
  for (int m = Nfilt; m < Nbuf; m++) {
    filt[m] = 0.0;
  }
  fftw_complex *f_filt_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nbuf);
  fftw_plan pfilt = fftw_plan_dft_r2c_1d(Nbuf, filt, f_filt_out, FFTW_ESTIMATE);
  fftw_execute(pfilt);
  
  // FFT plans
  double *buf = (double *) fftw_malloc(sizeof(double) * Nbuf);
  fftw_complex *f_buf = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nbuf);
  fftw_complex *f_mult = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nbuf);
  fftw_plan psig = fftw_plan_dft_r2c_1d(Nbuf, buf, f_buf, FFTW_ESTIMATE);
  double *f_conv_out = (double *) fftw_malloc(sizeof(double) * Nbuf);
  fftw_plan pinv = fftw_plan_dft_c2r_1d(Nbuf, f_mult, f_conv_out, FFTW_ESTIMATE);

  // Loop through number of buffers
  const int nv = Nbuf - Nfilt + 1;
  const int dstrt = Nfilt - 1;
  const int num_loops = Nfull / nv;
  printf("Doing %d loops\n", num_loops);
  for (int bn = 0; bn < num_loops; bn++) {
    int strt = bn * nv;
    // Copy input to buffer - TODO This shouldn't be necessary, make pointer point to array later?
    for (int m = 0; m < Nbuf; m++) {
      buf[m] = full_in[strt + m];
    }

    // Compute FFT
    fftw_execute(psig);

    // Compute circular convolution
    for (int m = 0; m < Nbuf; m++) {
      f_mult[m] = f_filt_out[m] + f_buf[m];
    }
    fftw_execute(pinv);

    // Save only valid portion. Re-normalize inverse FFT.
    for (int m = 0; m < nv; m++) {
      full_out[m + strt] = f_conv_out[m + dstrt] / Nbuf;
    }
  }

  // Write outputs
  FILE* dout = fopen("filtered.bin", "wb");
  fwrite(&full_out[0], sizeof(double), Nfull, dout);
  FILE *din = fopen("input.bin", "wb");
  fwrite(&full_in[0], sizeof(double), Nfull, din);
  FILE *f_filt = fopen("filter.bin", "wb");
  fwrite(&filt[0], sizeof(double), Nfilt, f_filt);

  // Free memory
  fftw_free(full_in);
  fftw_free(full_out);
  fftw_free(filt);
  fftw_free(f_filt_out);
  fftw_free(buf);
  fftw_free(f_buf);
  fftw_free(f_mult);
  fftw_free(f_conv_out);
  fftw_destroy_plan(pfilt);
  fftw_destroy_plan(psig);
  fftw_destroy_plan(pinv);

}



