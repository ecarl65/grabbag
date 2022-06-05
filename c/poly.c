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
  // Variables
  const int M = 4; // Downsample factor
  const int Nfull = 128 * M;  // Total number of samples
  const int Nds = Nfull / M;
  const int Nfilt = 16 * M + 1;
  const int Nfilt_half = (Nfilt - 1) >> 1;
  const double pi = acosf(-1);
  const double Fs = 10e3;
  const double Ts = 1.0 / Fs;
  const double fc = 100;
  const double ff = 2 * fc;  // Double the tone freq for the filter cutoff
  const int nrows = M;
  const int ncols = Nds;
  const int ncols_fft = ((Nfull / M) / 2 + 1);
  const int Nfft = ncols_fft * M;
  srand(time(NULL));

  // Arrays
  double *full_in = (double *) fftw_malloc(sizeof(double) * Nfull);
  double filt[Nfilt];
  double *filt_full = (double *) fftw_malloc(sizeof(double) * Nfull);
  double *ds_out = fftw_malloc(sizeof(double) * Nds);
  double *full_out = fftw_malloc(sizeof(double) * Nfull);
  fftw_complex *fft_in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nfft);
  fftw_complex *fft_filt = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nfft);
  double *conv_out = (double *) fftw_malloc(sizeof(double) * Nfull);
  fftw_complex *fft_mult = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nfft);

  // Make filter
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
  for (int m = 0; m < Nfilt; m++) filt[m] /= filt_sum;  // Normalize the filter

  // Make input sinusoid
  for (int m = 0; m < Nfull; m++) {
    full_in[m] = cos(2 * pi * fc * m * Ts) + coarse_gaussian() / 10.0;
  }

  // output
  for (size_t m = 0; m < Nfull; m++) full_out[m] = 0;
  for (size_t m = 0; m < Nds; m++) ds_out[m] = 0;

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
  n_size[0] = ncols_fft;
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
  for (int m = 0; m < Nfft; m++) {
    fft_mult[m] = fft_filt[m] * fft_in[m];
  }
  fftw_execute(pinv);

  // Perform summation and get downsampled output. Also normalize inverse fft.
  for (int m = 0; m < Nds; m++) {
    for (int rho = 0; rho < M; rho++) {
      ds_out[m] += full_out[m + rho * Nds] / Nfull;
    }
  }

  // Write outputs
  FILE* dout = fopen("filtered.bin", "wb");
  fwrite(&full_out[0], sizeof(double), Nfull, dout);
  fclose(dout);
  FILE *din = fopen("input.bin", "wb");
  fwrite(&full_in[0], sizeof(double), Nfull, din);
  fclose(din);
  FILE *f_filt = fopen("filter.bin", "wb");
  fwrite(&filt[0], sizeof(double), Nfilt, f_filt);
  fclose(f_filt);
  FILE* dsout = fopen("downsampled.bin", "wb");
  fwrite(&ds_out[0], sizeof(double), Nds, dsout);
  fclose(dsout);
  FILE* ddata = fopen("fftdata.bin", "wb");
  fwrite(&fft_in[0], sizeof(fftw_complex), Nfft, ddata);
  fclose(ddata);
  FILE* dfilt = fopen("fftfilt.bin", "wb");
  fwrite(&fft_filt[0], sizeof(fftw_complex), Nfft, dfilt);
  fclose(dfilt);

  // Free memory
  fftw_free(full_in);
  fftw_free(filt_full);
  fftw_free(ds_out);
  fftw_free(full_out);
  fftw_free(fft_in);
  fftw_free(fft_filt);
  fftw_free(conv_out);
  fftw_free(fft_mult);
  fftw_destroy_plan(pfilt);
  fftw_destroy_plan(psig);
  fftw_destroy_plan(pinv);
  fftw_cleanup();
}



