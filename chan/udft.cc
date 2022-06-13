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
#include <vector>
#include <fftw3.h>
#include <stdexcept>
#include "udft.hh"

// {{{ poly_filt_design
void UDFT::poly_filt_design()
{
  const int Nfilt_half = (n_filt - 1) >> 1;

  // Create plan
  fftwf_plan pfilt = fftwf_plan_many_dft_r2c(filt_c.rank, filt_c.n_size, filt_c.howmany, filt_full,
      filt_c.inembed, filt_c.istride, filt_c.idist, (fftwf_complex *) fft_filt,
      filt_c.onembed, filt_c.ostride, filt_c.odist, FFTW_MEASURE);

  // Make filter
  int outidx = 0;
  float filt_sum = 0;
  for (int m = -Nfilt_half; m <= Nfilt_half; m++) {
    // Sinc low pass filter at twice the SOI freq
    if (m == 0) {
      filt[outidx] = 2 * f_cutoff / samp_rate;
    } else {
      filt[outidx] = sin(2 * M_PI * f_cutoff * m / samp_rate) / (m * M_PI);
    }

    // Hamming window and normalization
    // filt[outidx] *= (0.54 - 0.46 * cos(2 * M_PI * outidx / n_filt));
    filt[outidx] *= ((25.0 / 46.0) - (21.0 / 46.0) * cos(2 * M_PI * outidx / n_filt));
    filt_sum += filt[outidx++];
  }
  for (int m = 0; m < n_filt; m++) filt[m] /= filt_sum;  // Normalize the filter

  // Filter polyphase decomposition
  for (int n = 0; n < n_cols; n++) {
    for (int rho = 0; rho < downsamp; rho++) {
      int inidx = n * downsamp - rho;
      int outidx = rho * n_cols + n;
      if (inidx < 0 || inidx >= n_filt) {
        filt_full[outidx] = 0;
      } else {
        filt_full[outidx] = filt[inidx];
      }
      outidx++;
    }
  }

  // Filter FFT
  fftwf_execute(pfilt);
  fftwf_destroy_plan(pfilt);
} // }}}

// {{{ UDFT
UDFT::UDFT(int downsamp, int n_full, int n_filt, float samp_rate, bool write, bool debug) : 
  downsamp(downsamp), n_full(n_full), n_filt(n_filt), samp_rate(samp_rate), write(write), debug(debug) 
{
  // Error checking on input
  if ((n_filt - 1) % (2 * downsamp) != 0) {
    throw std::invalid_argument("Filter length should be a multiple of 2x the downsample factor plus one\n");
  }
  if (downsamp % 2 != 0) {
    throw std::invalid_argument("Downsample should be an even amount\n");
  }

  // Buffer sizes
  n_buffer = (n_filt - 1) * 8;  // Size of the buffers, for now fix at 8x the filter size
  n_cols = n_buffer / downsamp;
  n_cols_fft = (n_cols / 2 + 1);
  n_rows_fft = downsamp / 2 + 1;
  n_fft_h = n_cols_fft * downsamp;
  n_fft_v = n_rows_fft * n_cols;
  n_out = n_full / downsamp * n_rows_fft;
  n_delay = (n_filt - 1) >> 1;
  n_delay_r = n_delay / downsamp;
  n_delay_samp = n_delay_r * n_rows_fft;
  idx_out_valid_r = (n_filt - 1) / downsamp;
  idx_out_valid_samp = idx_out_valid_r * n_rows_fft;
  n_in_valid = n_buffer - n_filt + 1;
  n_out_valid_r = n_in_valid / downsamp;
  n_out_valid_samp = n_out_valid_r * n_rows_fft;

  // Frequency and time constants
  samp_period = 1.0 / samp_rate;
  f_cutoff = samp_rate / (2 * downsamp);

  if (n_full < n_buffer) {
    throw std::invalid_argument("Input array length should be larger than 8x filter length\n");
  }

  if (debug) {
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
    printf("Starting output row of valid samples: %d\n", n_out_valid_r);
    printf("Number of valid output samples per buffer: %d\n", n_out_valid_samp);
    printf("Sample rate: %e\n", samp_rate);
    printf("Sample period: %e\n", samp_period);
    printf("Cutoff frequency: %e\n", f_cutoff);
  }

  // Forward FFT of the data
  fwd_c.n_size[0] = n_cols;
  fwd_c.rank = 1;
  fwd_c.howmany = downsamp;
  fwd_c.idist = 1;
  fwd_c.odist = n_cols_fft;
  fwd_c.istride = downsamp;
  fwd_c.ostride = 1;
  fwd_c.inembed = NULL;
  fwd_c.onembed = NULL;

  // Fwd FFT of the polyphase filter
  filt_c.n_size[0] = n_cols;
  filt_c.rank = 1;
  filt_c.howmany = downsamp;
  filt_c.idist = n_cols;
  filt_c.odist = n_cols_fft;
  filt_c.istride = 1;
  filt_c.ostride = 1;
  filt_c.inembed = NULL;
  filt_c.onembed = NULL;

  // Inverse FFT after multiplying the data and filter in the freq domain.
  inv_c.n_size[0] = n_cols;
  inv_c.rank = 1;
  inv_c.howmany = downsamp;
  inv_c.idist = n_cols_fft;
  inv_c.odist = n_cols;
  inv_c.istride = 1;
  inv_c.ostride = 1;
  inv_c.inembed = NULL;
  inv_c.onembed = NULL;

  // Transpose output FFT that's across channels and does the modulating.
  // The transpose makes it easier to skip the invalid samples and delayed samples.
  col_c.n_size[0] = downsamp;
  col_c.rank = 1;
  col_c.howmany = n_cols;
  col_c.idist = 1;
  col_c.odist = n_rows_fft;
  col_c.istride = n_cols;
  col_c.ostride = 1;
  col_c.inembed = NULL;
  col_c.onembed = NULL;

  // Arrays
  filt = fftwf_alloc_real(n_filt);
  buffer_in = fftwf_alloc_real(n_buffer);
  filt_full = fftwf_alloc_real(n_buffer);
  fft_in = reinterpret_cast<std::complex<float>*>(fftwf_alloc_complex(n_fft_h));
  fft_filt = reinterpret_cast<std::complex<float>*>(fftwf_alloc_complex(n_fft_h));
  conv_out = fftwf_alloc_real(n_buffer);
  fft_mult = reinterpret_cast<std::complex<float>*>(fftwf_alloc_complex(n_fft_h));
  udft = reinterpret_cast<std::complex<float>*>(fftwf_alloc_complex(n_fft_v));

  // Data polyphase FFT
  psig = fftwf_plan_many_dft_r2c(fwd_c.rank, fwd_c.n_size, fwd_c.howmany, buffer_in,
      fwd_c.inembed, fwd_c.istride, fwd_c.idist,
      reinterpret_cast<fftwf_complex*>(fft_in),
      fwd_c.onembed, fwd_c.ostride, fwd_c.odist, FFTW_MEASURE);

  // IFFT plan for convolution
  pinv = fftwf_plan_many_dft_c2r(inv_c.rank, inv_c.n_size, inv_c.howmany,
      reinterpret_cast<fftwf_complex*>(fft_mult),
      inv_c.inembed, inv_c.istride, inv_c.idist, conv_out,
      inv_c.onembed, inv_c.ostride, inv_c.odist, FFTW_MEASURE);

  // Perform FFT down columns to get channelized output. Output may be transposed.
  pudft = fftwf_plan_many_dft_r2c(col_c.rank, col_c.n_size, col_c.howmany, conv_out,
      col_c.inembed, col_c.istride, col_c.idist,
      reinterpret_cast<fftwf_complex*>(udft),
      col_c.onembed, col_c.ostride, col_c.odist, FFTW_MEASURE);

  // Design filter and do polyphase decomposition and FFT of filter
  poly_filt_design();

}
// }}}

// {{{ ~UDFT
UDFT::~UDFT() {
  // Free allocated arrays
  fftwf_free(filt);
  fftwf_free(buffer_in);
  fftwf_free(filt_full);
  fftwf_free(fft_in);
  fftwf_free(fft_filt);
  fftwf_free(conv_out);
  fftwf_free(fft_mult);
  fftwf_free(udft);

  // Do FFTW cleanup
  fftwf_destroy_plan(psig);
  fftwf_destroy_plan(pinv);
  fftwf_destroy_plan(pudft);
  fftw_cleanup_threads();
  fftwf_cleanup();
}
// }}}

// {{{ run
std::vector<std::vector<std::complex<float>>> UDFT::run(float *indata)
{
  // Initialize output
  int n_out_rows = n_full / downsamp;
  int n_out_cols = n_rows_fft;
  std::vector<std::vector<std::complex<float>>> full_out(n_out_rows, std::vector<std::complex<float>> (n_out_cols));
  for (int m = 0; m < n_delay_r; m++) {
    for (int n = 0; n < n_out_cols; n++) {
      full_out[m][n] = 0;
    }
  }

  // Move this to a function
  const int n_loops = (int) ceil((float) n_full / n_in_valid);
  if (debug) printf("Number of loops: %d\n", n_loops);
  for (int idx = 0; idx < n_loops; idx++) {
    int in_start = n_in_valid * idx;
    int out_start_r = n_in_valid / downsamp * idx;
    if (debug) printf("Input index range: [%d, %d)\n", in_start, in_start + n_cols * downsamp);

    // Forward FFT of this buffer of data. If last buffer do zero padding.
    if (n_full - in_start < n_buffer) {
      for (int m = 0; m < n_full - in_start; m++) buffer_in[m] = indata[in_start + m];
      for (int m = n_full - in_start; m < n_buffer; m++) buffer_in[m] = 0;
      fftwf_execute_dft_r2c(psig, buffer_in, reinterpret_cast<fftwf_complex*>(fft_in));
    } else {
      fftwf_execute_dft_r2c(psig, &indata[in_start], reinterpret_cast<fftwf_complex*>(fft_in));
    }

    // Compute circular convolution through multiplying in frequency domain.
    for (int m = 0; m < n_fft_h; m++) fft_mult[m] = fft_filt[m] * fft_in[m];

    // Perform inverse FFT to get real data out of convolution
    fftwf_execute(pinv);

    fftwf_execute(pudft);

    // Save only valid portion. Re-normalize inverse FFT.
    if (debug) printf("Copying from UDFT row %d to output row %d\n", idx_out_valid_r, out_start_r + n_delay_r);
    for (int r = 0; r < n_out_valid_r; r++) {
      for (int c = 0; c < n_out_cols; c++) {
        int m = n_out_cols * r + c;
        std::cout << "r = " << r << ", c = " << c << ", m = " << m << ", out_start_r = " << out_start_r 
          << ", n_delay_r = " << n_delay_r << ", idx_out_valid_samp = " << idx_out_valid_samp << 
          ", n_out_valid_r = " << n_out_valid_r << ", n_out_cols = " << n_out_cols << 
          ", out_row = " << r + out_start_r + n_delay_r << ", in_sample = " << m + idx_out_valid_samp << "\n";
        if (r + out_start_r + n_delay_r >= n_out_rows) break;
        full_out[r + out_start_r + n_delay_r][c] = udft[m + idx_out_valid_samp] / (float) n_cols;
      }
    }
  } // End loop

  // Write outputs
  if (write) {
    write_out("filtered.bin", (void *) &conv_out[0], sizeof(float), n_buffer);
    write_out("input.bin", (void *) &indata[0], sizeof(float), n_full);
    write_out("filter.bin", (void *) &filt[0], sizeof(float), n_filt);
    write_out("onebuffer.bin", (void *) &udft[0], sizeof(fftwf_complex), n_fft_v);
    // write_append("channelized.bin", (void *) &full_out[m], sizeof(fftwf_complex), n_out);
    for (size_t m = 0; m < full_out.size(); m++) {
      write_append("channelized.bin", (void *) &full_out[m][0], sizeof(fftwf_complex), full_out[m].size());
    }
    write_out("fftdata.bin", (void *) &fft_in[0], sizeof(fftwf_complex), n_fft_h);
    write_out("fftfilt.bin", (void *) &fft_filt[0], sizeof(fftwf_complex), n_fft_h);
  }

  return full_out;
}
// }}}

