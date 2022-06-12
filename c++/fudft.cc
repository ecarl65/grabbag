/*
 * =====================================================================================
 *
 *       Filename:  fudft.c
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
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <tgmath.h>
#include <fftw3.h>
#include <string.h>
#include "utils.h"

int channelizer(const int downsamp, const int n_full, const int n_filt, const float samp_rate) {
  // Error checking on input
  if ((n_filt - 1) % (2 * downsamp) != 0) {
    printf("Filter length should be a multiple of 2x the downsample factor plus one\n");
    return 0;
  }
  if (downsamp % 2 != 0) {
    printf("Downsample should be an even amount\n");
    return 0;
  }

  // Buffer sizes
  const int n_buffer = (n_filt - 1) * 8;  // Size of the buffers, for now fix at 8x the filter size
  const int n_cols = n_buffer / downsamp;
  const int n_cols_fft = (n_cols / 2 + 1);
  const int n_rows_fft = downsamp / 2 + 1;
  const int n_fft_h = n_cols_fft * downsamp;
  const int n_fft_v = n_rows_fft * n_cols;
  const int n_out = n_full / downsamp * n_rows_fft;
  const int n_delay = (n_filt - 1) >> 1;
  const int n_delay_r = n_delay / downsamp;
  const int n_delay_samp = n_delay_r * n_rows_fft;
  const int idx_out_valid_r = (n_filt - 1) / downsamp;
  const int idx_out_valid_samp = idx_out_valid_r * n_rows_fft;
  const int n_in_valid = n_buffer - n_filt + 1;
  const int n_out_valid = n_in_valid / downsamp * n_rows_fft;

  if (n_full < n_buffer) {
    printf("Input array length should be larger than 8x filter length\n");
    return 0;
  }

#ifdef DEBUG
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
  printf("Output row of valid overlap/save data: %d\n", idx_out_valid_r);
  printf("Output sample of valid overlap/save data: %d\n", idx_out_valid_samp);
  printf("Number of valid input samples per buffer: %d\n", n_in_valid);
  printf("Number of valid output samples per buffer: %d\n", n_out_valid);
#endif
  
  // Frequency and time constants
  const float samp_period = 1.0 / samp_rate;
  const float chirp_period = n_full * samp_period / 2;
  const float f_cutoff = samp_rate / (2 * downsamp);

#ifdef DEBUG
  printf("Sample rate: %e\n", samp_rate);
  printf("Sample period: %e\n", samp_period);
  printf("Chirp period: %e\n", chirp_period);
  printf("Cutoff frequency: %e\n", f_cutoff);
#endif

  // FFT parameters: n_size, rank, howmany, idist, odist, istride, ostride, inembed, onembed
  struct fft_config fwd_c = { {n_cols}, 1, downsamp, 1, n_cols_fft, downsamp, 1, NULL, NULL };
  struct fft_config filt_c = { {n_cols}, 1, downsamp, n_cols, n_cols_fft, 1, 1, NULL, NULL };
  struct fft_config inv_c = { {n_cols}, 1, downsamp, n_cols_fft, n_cols, 1, 1, NULL, NULL };
  struct fft_config col_c = {
    {downsamp}, 1, n_cols, 1, n_rows_fft, n_cols, 1, NULL, NULL   // Transposed version - Makes overlap/save easier
  };

  // Arrays
  float *full_in = fftwf_alloc_real(n_full);
  std::complex<float> *full_out = reinterpret_cast<std::complex<float>*>(fftwf_alloc_complex(n_out));
  float filt[n_filt];  // Short, non-zero taps only, just leave on the stack. Will copy to below before FFT.
  float *buffer_in = fftwf_alloc_real(n_buffer);
  float *filt_full = fftwf_alloc_real(n_buffer);
  std::complex<float> *fft_in = reinterpret_cast<std::complex<float>*>(fftwf_alloc_complex(n_fft_h));
  std::complex<float> *fft_filt = reinterpret_cast<std::complex<float>*>(fftwf_alloc_complex(n_fft_h));
  float *conv_out = fftwf_alloc_real(n_buffer);
  std::complex<float> *fft_mult = reinterpret_cast<std::complex<float>*>(fftwf_alloc_complex(n_fft_h));
  std::complex<float> *udft = reinterpret_cast<std::complex<float>*>(fftwf_alloc_complex(n_fft_v));

  for (size_t m = 0; m < n_delay_samp; m++) full_out[m] = 0;

  // Make input chirp
  make_chirp(full_in, n_full, samp_rate, chirp_period);

  // Design filter and do polyphase decomposition and FFT of filter
  poly_filt_design(n_filt, f_cutoff, samp_rate, &filt[0], filt_full, n_cols, downsamp, filt_c, fft_filt);

  // Data polyphase FFT
  fftwf_plan psig = fftwf_plan_many_dft_r2c(fwd_c.rank, fwd_c.n_size, fwd_c.howmany, buffer_in,
                                          fwd_c.inembed, fwd_c.istride, fwd_c.idist,
                                          reinterpret_cast<fftwf_complex*>(fft_in),
                                          fwd_c.onembed, fwd_c.ostride, fwd_c.odist, FFTW_ESTIMATE);

  // IFFT plan for convolution
  fftwf_plan pinv = fftwf_plan_many_dft_c2r(inv_c.rank, inv_c.n_size, inv_c.howmany,
                                          reinterpret_cast<fftwf_complex*>(fft_mult),
                                          inv_c.inembed, inv_c.istride, inv_c.idist, conv_out,
                                          inv_c.onembed, inv_c.ostride, inv_c.odist, FFTW_ESTIMATE);

  // Perform FFT down columns to get channelized output. Output may be transposed.
  fftwf_plan pudft = fftwf_plan_many_dft_r2c(col_c.rank, col_c.n_size, col_c.howmany, conv_out,
                                           col_c.inembed, col_c.istride, col_c.idist,
                                           reinterpret_cast<fftwf_complex*>(udft),
                                           col_c.onembed, col_c.ostride, col_c.odist, FFTW_ESTIMATE);

  // Move this to a function
  const int n_loops = (int) ceil((float) n_full / n_in_valid);
#ifdef DEBUG
  printf("Number of loops: %d\n", n_loops);
#endif
  for (size_t idx = 0; idx < n_loops; idx++) {  // TODO Get right number of buffers
    int in_start = n_in_valid * idx;
    int out_start = n_out_valid * idx;
#ifdef DEBUG
    printf("Input index range: [%d, %d)\n", in_start, in_start + n_cols * downsamp);
#endif

    // Forward FFT of this buffer of data. If last buffer do zero padding.
    if (n_full - in_start < n_buffer) {
      for (size_t m = 0; m < n_full - in_start; m++) buffer_in[m] = full_in[in_start + m];
      for (size_t m = n_full - in_start; m < n_buffer; m++) buffer_in[m] = 0;
      fftwf_execute_dft_r2c(psig, buffer_in, reinterpret_cast<fftwf_complex*>(fft_in));
    } else {
      fftwf_execute_dft_r2c(psig, &full_in[in_start], reinterpret_cast<fftwf_complex*>(fft_in));
    }

    // Compute circular convolution through multiplying in frequency domain.
    for (int m = 0; m < n_fft_h; m++) fft_mult[m] = fft_filt[m] * fft_in[m];

    // Perform inverse FFT to get real data out of convolution
    fftwf_execute(pinv);

    fftwf_execute(pudft);

    // Save only valid portion. Re-normalize inverse FFT.
#ifdef DEBUG
    printf("Copying from UDFT %d to output %d\n", idx_out_valid_samp, out_start + n_delay_samp);
#endif
    for (int m = 0; m < n_out_valid; m++) {
      if (m + out_start + n_delay_samp >= n_out) break;
      full_out[m + out_start + n_delay_samp] = udft[m + idx_out_valid_samp] / (float) n_cols;
    }

  } // End loop

  // Write outputs
  write_out("filtered.bin", (void *) &conv_out[0], sizeof(float), n_buffer);
  write_out("input.bin", (void *) &full_in[0], sizeof(float), n_full);
  write_out("filter.bin", (void *) &filt[0], sizeof(float), n_filt);
  write_out("onebuffer.bin", (void *) &udft[0], sizeof(fftwf_complex), n_fft_v);
  write_out("channelized.bin", (void *) &full_out[0], sizeof(fftwf_complex), n_out);
  write_out("fftdata.bin", (void *) &fft_in[0], sizeof(fftwf_complex), n_fft_h);
  write_out("fftfilt.bin", (void *) &fft_filt[0], sizeof(fftwf_complex), n_fft_h);

  // Free memory
  fftwf_free(full_in);
  fftwf_free(full_out);
  fftwf_free(filt_full);
  fftwf_free(buffer_in);
  fftwf_free(udft);
  fftwf_free(fft_in);
  fftwf_free(fft_filt);
  fftwf_free(conv_out);
  fftwf_free(fft_mult);
  fftwf_destroy_plan(psig);
  fftwf_destroy_plan(pinv);
  fftwf_destroy_plan(pudft);
  fftw_cleanup_threads();
  fftwf_cleanup();

  return 1;

}

int main(int argc, char **argv) {
  // Values on which all others depend
  int downsamp = 8; // Downsample factor
  int filt_ord = 8;
  int full_ord = 256;
  float samp_rate = 10e3;

  int opt;

  // put ':' in the starting of the string so that program can distinguish between '?' and ':' 
  while((opt = getopt(argc, argv, ":d:f::n::s::h::")) != -1) 
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
        printf("\th - This help\n");
        exit(EXIT_FAILURE);
        break;
      case 'd': 
        downsamp = atoi(optarg);
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

  int valid = channelizer(downsamp, n_full, n_filt, samp_rate);

  if (!valid) {
    printf("Processing FAILED\n");
    exit(EXIT_FAILURE);
  }
}



