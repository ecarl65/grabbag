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

double coarse_gaussian() {
  double out = 0;
  for (int m = 0; m < 12; m++) {
    out += (double) (rand() % 100);
  }
  out /= 100.0;
  out -= 6.0;
}

int main(int argc, char **argv) {
  // Make input sinusoid
  const int Nbuf = 1024;
  const int Nbufs = 20;
  const int Nfull = Nbuf * Nbufs;
  const double pi = acosf(-1);
  const double Fs = 10e3;
  const double Ts = 1.0 / Fs;
  const double fc = 100;
  const double ff = 2 * fc;
  double *full_in = (double *) fftw_malloc(sizeof(double) * Nfull);
  srand(time(NULL));
  for (int m = 0; m < Nfull; m++) {
    full_in[m] = cos(2 * pi * fc * m * Ts) + coarse_gaussian() / 10.0;
  }
  double *full_out = fftw_malloc(sizeof(double) * Nfull);
  memset(&full_out[0], 0, sizeof(full_out[0]) * Nfull);

  // Make filter
  const int Nfilt = 129;
  const int Nfilt_half = (Nfilt - 1) >> 1;
  double *filt = (double *) fftw_malloc(sizeof(double) * Nbuf);
  int outidx = 0;
  double filt_sum = 0;
  for (int m = -Nfilt_half; m <= Nfilt_half; m++) {
    // Sinc low pass filter at twice the SOI freq
    if (m == 0) {
      filt[outidx] = 2 * ff / Fs;
    } else {
      filt[outidx] = sin(2 * pi * ff * m / Fs) / (m * pi);
    }

    // Hamming window and normalization
    filt[outidx] *= (0.54 - 0.46 * cos(2 * pi * outidx / Nfilt));
    filt_sum += filt[outidx++];
  }
  for (int m = 0; m < Nbuf; m++) {
    if (m < Nfilt) {
      filt[m] /= filt_sum;
    } else {
      filt[m] = 0.0f;
    }
  }
  fftw_complex *f_filt_out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nbuf);
  fftw_plan pfilt = fftw_plan_dft_r2c_1d(Nbuf, filt, f_filt_out, FFTW_ESTIMATE);
  fftw_execute(pfilt);
  
  // FFT plans
  double *buf = (double *) fftw_malloc(sizeof(double) * Nbuf);
  fftw_complex *f_buf = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nbuf);
  fftw_plan psig = fftw_plan_dft_r2c_1d(Nbuf, buf, f_buf, FFTW_ESTIMATE);
  double *conv_out = (double *) fftw_malloc(sizeof(double) * Nbuf);
  fftw_plan pinv = fftw_plan_dft_c2r_1d(Nbuf, f_buf, conv_out, FFTW_ESTIMATE);

  // TODO Use wisdom

  // Loop through number of buffers
  const int nv = Nbuf - Nfilt + 1;
  const int dstrt = Nfilt - 1;
  const int num_loops = Nfull / nv;
  const int delay = dstrt >> 1;
  const int ds = 4;
  const int Nds = Nfull / ds;
  double *ds_out = fftw_malloc(sizeof(double) * Nds);
  printf("Doing %d loops\n", num_loops);
  for (int bn = 0; bn < num_loops; bn++) {
    // Forward FFT
    int strt = bn * nv;
    fftw_execute_dft_r2c(psig, &full_in[strt], f_buf);

    // Compute circular convolution
    for (int m = 0; m < Nbuf; m++) {
      f_buf[m] *= f_filt_out[m];
    }
    fftw_execute(pinv);

    // Save only valid portion. Re-normalize inverse FFT.
    for (int m = 0; m < nv; m++) {
      int outidx = m + strt + delay;
      full_out[outidx] = conv_out[m + dstrt] / Nbuf;
      if (outidx % ds == 0) {
        ds_out[outidx / ds] = full_out[outidx];
      }
    }
  }

  // Perform downsampling
  for (int m = 0; m < Nfull; m++) {
    if (m % ds == 0) {
      /* ds_out[m/ds] = full_out[m]; */
    }
  }

  // Write outputs
  FILE* dsout = fopen("downsampled.bin", "wb");
  fwrite(&ds_out[0], sizeof(double), Nds, dsout);
  fclose(dsout);
  FILE* dout = fopen("filtered.bin", "wb");
  fwrite(&full_out[0], sizeof(double), Nfull, dout);
  fclose(dout);
  FILE *din = fopen("input.bin", "wb");
  fwrite(&full_in[0], sizeof(double), Nfull, din);
  fclose(din);
  FILE *f_filt = fopen("filter.bin", "wb");
  fwrite(&filt[0], sizeof(double), Nfilt, f_filt);
  fclose(f_filt);

  // Free memory
  fftw_free(full_in);
  fftw_free(full_out);
  fftw_free(ds_out);
  fftw_free(filt);
  fftw_free(f_filt_out);
  fftw_free(buf);
  fftw_free(f_buf);
  fftw_free(conv_out);
  fftw_destroy_plan(pfilt);
  fftw_destroy_plan(psig);
  fftw_destroy_plan(pinv);
  fftw_cleanup();
}



