/*
 * =====================================================================================
 *
 *       Filename:  fudft.hh
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

// {{{ poly_filt_design
template <class T>
void poly_filt_design(int Nfilt, T fcutoff, T Fs, T* filt, T* filt_full, int ncols, int downsamp,
                      const struct fft_config filt_cfg, std::complex<T> *fft_filt) 
{
  T Ts = 1 / Fs;
  const int Nfilt_half = (Nfilt - 1) >> 1;

  // Make filter
  int outidx = 0;
  T filt_sum = 0;
  for (int m = -Nfilt_half; m <= Nfilt_half; m++) {
    // Sinc low pass filter at twice the SOI freq
    if (m == 0) {
      filt[outidx] = 2 * fcutoff / Fs;
    } else {
      filt[outidx] = sin(2 * M_PI * fcutoff * m / Fs) / (m * M_PI);
    }

    // Hamming window and normalization
    // filt[outidx] *= (0.54 - 0.46 * cos(2 * M_PI * outidx / Nfilt));
    filt[outidx] *= ((25.0 / 46.0) - (21.0 / 46.0) * cos(2 * M_PI * outidx / Nfilt));
    filt_sum += filt[outidx++];
  }
  for (int m = 0; m < Nfilt; m++) filt[m] /= filt_sum;  // Normalize the filter

  // Filter polyphase decomposition
  for (int n = 0; n < ncols; n++) {
    for (int rho = 0; rho < downsamp; rho++) {
      int inidx = n * downsamp - rho;
      int outidx = rho * ncols + n;
      if (inidx < 0 || inidx >= Nfilt) {
        filt_full[outidx] = 0;
      } else {
        filt_full[outidx] = filt[inidx];
      }
      outidx++;
    }
  }

  // Filter FFT
  fftwf_plan pfilt = fftwf_plan_many_dft_r2c(filt_cfg.rank, filt_cfg.n_size, filt_cfg.howmany, filt_full,
                                           filt_cfg.inembed, filt_cfg.istride, filt_cfg.idist, (fftwf_complex *) fft_filt,
                                           filt_cfg.onembed, filt_cfg.ostride, filt_cfg.odist, FFTW_ESTIMATE);
  fftwf_execute(pfilt);
  fftwf_destroy_plan(pfilt);
}
// }}} 
