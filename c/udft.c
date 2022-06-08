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
  const int downsamp = 8; // Downsample factor

  // Buffer sizes
  const int n_full = 256 * downsamp;  // Total number of samples
  const int n_filt = 8 * downsamp + 1;
  const int n_buffer = (n_filt - 1) * 8;  // Size of the buffers
  const int n_cols = n_buffer / downsamp;
  const int n_cols_fft = ((n_buffer / downsamp) / 2 + 1);
  const int n_rows_fft = downsamp / 2 + 1;
  const int n_fft_h = n_cols_fft * downsamp;
  const int n_fft_v = n_rows_fft * n_cols;
  const int n_out = n_full / downsamp * n_rows_fft;
  const int n_delay = (n_filt - 1) >> 1;
  const int n_delay_r = n_delay / downsamp;
  const int n_delay_samp = n_delay_r * n_rows_fft;
  const int out_valid_r = (n_filt - 1) / downsamp;
  const int out_valid_samp = out_valid_r * n_rows_fft;
  const int n_valid = (n_buffer - n_filt + 1) / downsamp * n_rows_fft;
  /* const int n_channel_delay = n_delay / downsamp * n_rows_fft; */

  printf("Downsample amount: %d\n", downsamp);
  printf("Number of input samples: %d\n", n_full);
  printf("Filter length: %d\n", n_filt);
  printf("Size of buffer: %d\n", n_buffer);
  printf("Number of samples per output channel per buffer: %d\n", n_cols);
  printf("Number of samples per output channel per buffer in freq domain: %d\n", n_cols_fft);
  printf("Number of output channels: %d\n", n_rows_fft);
  printf("Number of samples in FFT of data and filter: %d\n", n_fft_h);
  printf("Number of samples in modulating FFT across channels: %d\n", n_fft_v);
  printf("Total number of output samples: %d\n", n_out);
  printf("Filter delay in samples at input rate: %d\n", n_delay);
  printf("Output row for zero group delay of filter: %d\n", n_delay_r);
  printf("Output sample for zero group delay of filter: %d\n", n_delay_samp);
  printf("Output row of valid overlap/save data: %d\n", out_valid_r);
  printf("Output sample of valid overlap/save data: %d\n", out_valid_samp);
  printf("Number of valid output samples: %d\n", n_valid);
  
  // Frequency and time constants
  const double samp_rate = 10e3;
  const double samp_period = 1.0 / samp_rate;
  const double chirp_period = n_full * samp_period / 2;
  const double f_cutoff = samp_rate / (2 * downsamp);

  printf("Sample rate: %e\n", samp_rate);
  printf("Sample period: %e\n", samp_period);
  printf("Chirp period: %e\n", chirp_period);
  printf("Cutoff frequency: %e\n", f_cutoff);

  // FFT parameters: n_size, rank, howmany, idist, odist, istride, ostride, inembed, onembed
  struct fft_config fwd_c = {
    {n_cols}, 1, downsamp, 1, n_cols_fft, downsamp, 1, NULL, NULL
  };
  struct fft_config filt_c = {
    {n_cols}, 1, downsamp, n_cols, n_cols_fft, 1, 1, NULL, NULL
  };
  struct fft_config inv_c = {
    {n_cols}, 1, downsamp, n_cols_fft, n_cols, 1, 1, NULL, NULL
  };
  struct fft_config col_c = {
    // {downsamp}, 1, n_cols, 1, 1, n_cols, n_cols, NULL, NULL       // Non-transposed
    {downsamp}, 1, n_cols, 1, n_rows_fft, n_cols, 1, NULL, NULL   // Transposed version - Makes overlap/save easier
  };

  // Arrays
  double *full_in = fftw_alloc_real(n_full);
  fftw_complex *full_out = fftw_alloc_complex(n_out);
  double filt[n_filt];  // Short, non-zero taps only, just leave on the stack. Will copy to below before FFT.
  double *buffer_in = fftw_alloc_real(n_buffer);
  double *filt_full = fftw_alloc_real(n_buffer);
  fftw_complex *fft_in = fftw_alloc_complex(n_fft_h);
  fftw_complex *fft_filt = fftw_alloc_complex(n_fft_h);
  double *conv_out = fftw_alloc_real(n_buffer);
  fftw_complex *fft_mult = fftw_alloc_complex(n_fft_h);
  fftw_complex *udft = fftw_alloc_complex(n_fft_v);

  // Design filter and do polyphase decomposition
  poly_filt_design(n_filt, f_cutoff, samp_rate, &filt[0], filt_full, n_cols, downsamp);

  // Make input chirp
  make_chirp(full_in, n_full, samp_rate, chirp_period);

  // Data polyphase FFT
  fftw_plan psig = fftw_plan_many_dft_r2c(fwd_c.rank, fwd_c.n_size, fwd_c.howmany, buffer_in,
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

  // Perform FFT down columns to get channelized output. Output may be transposed.
  fftw_plan pudft = fftw_plan_many_dft_r2c(col_c.rank, col_c.n_size, col_c.howmany, conv_out,
                                           col_c.inembed, col_c.istride, col_c.idist, udft,
                                           col_c.onembed, col_c.ostride, col_c.odist, FFTW_ESTIMATE);

  // TODO: this section will be in the loop
  const int num_loops = (n_full - n_cols) / n_valid;
  printf("Number of loops: %d\n", num_loops);
  for (size_t idx = 0; idx < 4; idx++) {  // TODO Get right number of buffers
    int in_start = n_valid * idx;
    /* int in_start = n_fft_v * idx; */

    // Forward FFT of this buffer of data
    fftw_execute_dft_r2c(psig, &full_in[in_start], fft_in);

    // Compute circular convolution through multiplying in frequency domain.
    for (int m = 0; m < n_fft_h; m++) fft_mult[m] = fft_filt[m] * fft_in[m];

    // Perform inverse FFT to get real data out of convolution
    fftw_execute(pinv);

    fftw_execute(pudft);

    // TODO Copy udft to output, or append to output file
    // Save only valid portion. Re-normalize inverse FFT.
    /* for (int m = 0; m < n_fft_v; m++) { */
    printf("Copying from UDFT output to %d\n", in_start + n_delay_samp);
    for (int m = 0; m < n_valid; m++) {
      full_out[m + in_start + n_delay_samp] = udft[m + out_valid_samp] * 4;
      /* full_out[m + in_start] = udft[m] * 4;  // Why 4 ? Downsamp/2? */
    }

  } // End loop

  // Write outputs
  write_out("filtered.bin", (void *) &conv_out[0], sizeof(double), n_buffer);
  write_out("input.bin", (void *) &full_in[0], sizeof(double), n_full);
  write_out("filter.bin", (void *) &filt[0], sizeof(double), n_filt);
  write_out("onebuffer.bin", (void *) &udft[0], sizeof(fftw_complex), n_fft_v);
  write_out("channelized.bin", (void *) &full_out[0], sizeof(fftw_complex), n_out);
  write_out("fftdata.bin", (void *) &fft_in[0], sizeof(fftw_complex), n_fft_h);
  write_out("fftfilt.bin", (void *) &fft_filt[0], sizeof(fftw_complex), n_fft_h);

  // Free memory
  fftw_free(full_in);
  fftw_free(full_out);
  fftw_free(filt_full);
  fftw_free(buffer_in);
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



