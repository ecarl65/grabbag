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

template <class T>
class UDFT {
  public:
    // Variables
    int downsamp, n_full, n_filt;
    T samp_rate;
    int n_buffer, n_cols, n_cols_fft, n_rows_fft, n_fft_h, n_fft_v, n_out, n_delay, n_delay_r, n_delay_samp;
    int idx_out_valid_r, idx_out_valid_samp, n_in_valid, n_out_valid;
    T samp_period, chirp_period, f_cutoff;
    struct fft_config fwd_c, filt_c, inv_c, col_c;
    T *full_in, *filt, *buffer_in, *filt_full, *conv_out;
    std::complex<T> *full_out, *fft_in, *fft_filt, *fft_mult, *udft;
    fftwf_plan psig, pinv, pudft;
    bool debug, write;

    // {{{ poly_filt_design
    inline void poly_filt_design()
    {
      T Ts = 1 / samp_rate;
      const int Nfilt_half = (n_filt - 1) >> 1;

      // Make filter
      int outidx = 0;
      T filt_sum = 0;
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
      fftwf_plan pfilt = fftwf_plan_many_dft_r2c(filt_c.rank, filt_c.n_size, filt_c.howmany, filt_full,
          filt_c.inembed, filt_c.istride, filt_c.idist, (fftwf_complex *) fft_filt,
          filt_c.onembed, filt_c.ostride, filt_c.odist, FFTW_MEASURE);
      fftwf_execute(pfilt);
      fftwf_destroy_plan(pfilt);
    } // }}}

    // {{{ UDFT
    UDFT(int downsamp, int n_full, int n_filt, T samp_rate, bool write = false, bool debug = true) : 
      downsamp(downsamp), n_full(n_full), n_filt(n_filt), samp_rate(samp_rate), write(write), debug(debug) {
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
        n_out_valid = n_in_valid / downsamp * n_rows_fft;

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
          printf("Number of valid output samples per buffer: %d\n", n_out_valid);
        }

        // Frequency and time constants
        samp_period = 1.0 / samp_rate;
        chirp_period = n_full * samp_period / 2;
        f_cutoff = samp_rate / (2 * downsamp);

        if (debug) {
          printf("Sample rate: %e\n", samp_rate);
          printf("Sample period: %e\n", samp_period);
          printf("Chirp period: %e\n", chirp_period);
          printf("Cutoff frequency: %e\n", f_cutoff);
        }

        // FFT parameters: n_size, rank, howmany, idist, odist, istride, ostride, inembed, onembed
        fwd_c.n_size[0] = n_cols;
        fwd_c.rank = 1;
        fwd_c.howmany = downsamp;
        fwd_c.idist = 1;
        fwd_c.odist = n_cols_fft;
        fwd_c.istride = downsamp;
        fwd_c.ostride = 1;
        fwd_c.inembed = NULL;
        fwd_c.onembed = NULL;

        filt_c.n_size[0] = n_cols;
        filt_c.rank = 1;
        filt_c.howmany = downsamp;
        filt_c.idist = n_cols;
        filt_c.odist = n_cols_fft;
        filt_c.istride = 1;
        filt_c.ostride = 1;
        filt_c.inembed = NULL;
        filt_c.onembed = NULL;

        inv_c.n_size[0] = n_cols;
        inv_c.rank = 1;
        inv_c.howmany = downsamp;
        inv_c.idist = n_cols_fft;
        inv_c.odist = n_cols;
        inv_c.istride = 1;
        inv_c.ostride = 1;
        inv_c.inembed = NULL;
        inv_c.onembed = NULL;

        // Transpose
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
        full_in = fftwf_alloc_real(n_full);
        full_out = reinterpret_cast<std::complex<T>*>(fftwf_alloc_complex(n_out));
        filt = fftwf_alloc_real(n_filt);
        buffer_in = fftwf_alloc_real(n_buffer);
        filt_full = fftwf_alloc_real(n_buffer);
        fft_in = reinterpret_cast<std::complex<T>*>(fftwf_alloc_complex(n_fft_h));
        fft_filt = reinterpret_cast<std::complex<T>*>(fftwf_alloc_complex(n_fft_h));
        conv_out = fftwf_alloc_real(n_buffer);
        fft_mult = reinterpret_cast<std::complex<T>*>(fftwf_alloc_complex(n_fft_h));
        udft = reinterpret_cast<std::complex<T>*>(fftwf_alloc_complex(n_fft_v));

        for (size_t m = 0; m < n_delay_samp; m++) full_out[m] = 0;

        // Make input chirp
        make_chirp(full_in, n_full, samp_rate, chirp_period);

        // Design filter and do polyphase decomposition and FFT of filter
        poly_filt_design();

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
      }
    // }}}

    // {{{ ~UDFT
    ~UDFT() {
      // Free memory
      fftwf_free(full_in);
      fftwf_free(full_out);
      fftwf_free(filt_full);
      fftwf_free(buffer_in);
      fftwf_free(udft);
      fftwf_free(fft_in);
      fftwf_free(filt);
      fftwf_free(fft_filt);
      fftwf_free(conv_out);
      fftwf_free(fft_mult);
      fftwf_destroy_plan(psig);
      fftwf_destroy_plan(pinv);
      fftwf_destroy_plan(pudft);
      fftw_cleanup_threads();
      fftwf_cleanup();
    }
    // }}}

    // {{{ run
    inline void run()
    {
      // Move this to a function
      const int n_loops = (int) ceil((T) n_full / n_in_valid);
      if (debug) printf("Number of loops: %d\n", n_loops);
      for (size_t idx = 0; idx < n_loops; idx++) {
        int in_start = n_in_valid * idx;
        int out_start = n_out_valid * idx;
        if (debug) printf("Input index range: [%d, %d)\n", in_start, in_start + n_cols * downsamp);

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
        if (debug) printf("Copying from UDFT %d to output %d\n", idx_out_valid_samp, out_start + n_delay_samp);
        for (int m = 0; m < n_out_valid; m++) {
          if (m + out_start + n_delay_samp >= n_out) break;
          full_out[m + out_start + n_delay_samp] = udft[m + idx_out_valid_samp] / (T) n_cols;
        }

      } // End loop

      // Write outputs
      write_out("filtered.bin", (void *) &conv_out[0], sizeof(T), n_buffer);
      write_out("input.bin", (void *) &full_in[0], sizeof(T), n_full);
      write_out("filter.bin", (void *) &filt[0], sizeof(T), n_filt);
      write_out("onebuffer.bin", (void *) &udft[0], sizeof(fftwf_complex), n_fft_v);
      write_out("channelized.bin", (void *) &full_out[0], sizeof(fftwf_complex), n_out);
      write_out("fftdata.bin", (void *) &fft_in[0], sizeof(fftwf_complex), n_fft_h);
      write_out("fftfilt.bin", (void *) &fft_filt[0], sizeof(fftwf_complex), n_fft_h);
    }
    // }}}

};

