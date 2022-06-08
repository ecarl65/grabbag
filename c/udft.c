/*
 * =====================================================================================
 *
 *       Filename:  udft.c
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
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <tgmath.h>
#include <fftw3.h>
#include <string.h>
#include "utils.h"

int main(int argc, char **argv) {
  // Variables
  const int downsamp = 8; // Downsample factor
  const int n_full = 256 * downsamp;  // Total number of samples
  const int n_filt = 8 * downsamp + 1;
  const int n_buffer = (n_filt - 1) * 4;
  const int n_delay = (n_filt - 1) >> 1;
  const double samp_rate = 10e3;
  const double samp_period = 1.0 / samp_rate;
  const double chirp_period = n_full * samp_period / 2;
  const double f_cutoff = samp_rate / (2 * downsamp);
  const int ncols = n_full / downsamp;
  const int ncols_fft = ((n_full / downsamp) / 2 + 1);
  const int nrows_fft = downsamp / 2 + 1;
  const int n_channel_delay = n_delay / downsamp * nrows_fft;
  const int n_fft_h = ncols_fft * downsamp;
  const int n_fft_v = nrows_fft * ncols;

  // FFT parameters: n_size, rank, howmany, idist, odist, istride, ostride, inembed, onembed
  struct fft_config fwd_c = {
    {ncols}, 1, downsamp, 1, ncols_fft, downsamp, 1, NULL, NULL
  };
  struct fft_config filt_c = {
    {ncols}, 1, downsamp, ncols, ncols_fft, 1, 1, NULL, NULL
  };
  struct fft_config inv_c = {
    {ncols}, 1, downsamp, ncols_fft, ncols, 1, 1, NULL, NULL
  };
  struct fft_config col_c = {
    // {downsamp}, 1, ncols, 1, 1, ncols, ncols, NULL, NULL       // Non-transposed
    {downsamp}, 1, ncols, 1, nrows_fft, ncols, 1, NULL, NULL   // Transposed version - Makes O/S easier
  };

  // Arrays
  double *full_in = (double *) fftw_malloc(sizeof(double) * n_full);
  double filt[n_filt];  // Short, non-zero taps only, just leave on the stack. Will copy to below before FFT.
  double *filt_full = (double *) fftw_malloc(sizeof(double) * n_full);
  fftw_complex *fft_in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft_h);
  fftw_complex *fft_filt = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft_h);
  double *conv_out = (double *) fftw_malloc(sizeof(double) * n_full);
  fftw_complex *fft_mult = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft_h);
  fftw_complex *udft = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n_fft_v);

  // Design filter and do polyphase decomposition
  poly_filt_design(n_filt, f_cutoff, samp_rate, &filt[0], filt_full, ncols, downsamp);

  // Make input chirp
  make_chirp(full_in, n_full, samp_rate, chirp_period);

  // Data polyphase FFT
  fftw_plan psig = fftw_plan_many_dft_r2c(fwd_c.rank, fwd_c.n_size, fwd_c.howmany, full_in,
                                          fwd_c.inembed, fwd_c.istride, fwd_c.idist, fft_in,
                                          fwd_c.onembed, fwd_c.ostride, fwd_c.odist, FFTW_ESTIMATE);

  // Filter FFT
  fftw_plan pfilt = fftw_plan_many_dft_r2c(filt_c.rank, filt_c.n_size, filt_c.howmany, filt_full,
                                           filt_c.inembed, filt_c.istride, filt_c.idist, fft_filt,
                                           filt_c.onembed, filt_c.ostride, filt_c.odist, FFTW_ESTIMATE);
  fftw_execute(pfilt);
  
  // IFFT plan for convolution
  fftw_plan pinv = fftw_plan_many_dft_c2r(inv_c.rank, inv_c.n_size, inv_c.howmany, fft_mult,
                                          inv_c.inembed, inv_c.istride, inv_c.idist, conv_out,
                                          inv_c.onembed, inv_c.ostride, inv_c.odist, FFTW_ESTIMATE);

  // TODO: this section will be in the loop

  // Forward FFT of this buffer of data
  fftw_execute(psig);
  // fftw_execute_dft_r2c(psig, &full_in[idx], fft_in);

  // Compute circular convolution through multiplying in frequency domain.
  for (int m = 0; m < n_fft_h; m++) fft_mult[m] = fft_filt[m] * fft_in[m];

  // Perform inverse FFT to get real data out of convolution
  fftw_execute(pinv);

  // Perform FFT down columns to get channelized output. Output may be transposed.
  fftw_plan pudft = fftw_plan_many_dft_r2c(col_c.rank, col_c.n_size, col_c.howmany, conv_out,
                                           col_c.inembed, col_c.istride, col_c.idist, udft,
                                           col_c.onembed, col_c.ostride, col_c.odist, FFTW_ESTIMATE);
  fftw_execute(pudft);

  // TODO Copy udft to output, or append to output file

  // End loop

  // Write outputs
  write_out("filtered.bin", (void *) &conv_out[0], sizeof(double), n_full);
  write_out("input.bin", (void *) &full_in[0], sizeof(double), n_full);
  write_out("filter.bin", (void *) &filt[0], sizeof(double), n_filt);
  write_out("channelized.bin", (void *) &udft[n_channel_delay], sizeof(fftw_complex), n_fft_v - n_channel_delay);
  write_out("fftdata.bin", (void *) &fft_in[0], sizeof(fftw_complex), n_fft_h);
  write_out("fftfilt.bin", (void *) &fft_filt[0], sizeof(fftw_complex), n_fft_h);

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



