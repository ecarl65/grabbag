#include <iostream>
#include <complex>
#include <string>
#include <stdlib.h>
#include <stdint.h>
#include <tgmath.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

#define M_PI 3.14159265358979323846

// {{{ fft_config
struct fft_config {
  int n_size[1];
  int rank;
  int howmany;
  int idist;
  int odist;
  int istride;
  int ostride;
  int *inembed;
  int *onembed;
};
// }}}

// {{{ write_out
void write_out(std::string file, void *addr, size_t size, size_t numels) {
  FILE* fout = fopen(file.c_str(), "wb");
  fwrite(addr, size, numels, fout);
  fclose(fout);
}
// }}}

// {{{ write_append
void write_append(std::string file, void *addr, size_t size, size_t numels) {
  FILE* fout = fopen(file.c_str(), "ab");
  fwrite(addr, size, numels, fout);
  fclose(fout);
}
// }}}

// {{{ poly_filt_design_f
template <class T>
void poly_filt_design(int Nfilt, T fcutoff, T Fs, T* filt, T* filt_full, int ncols, int downsamp,
                      const struct fft_config filt_cfg, std::complex<T> *fft_filt) 
{
  float Ts = 1 / Fs;
  const int Nfilt_half = (Nfilt - 1) >> 1;

  // Make filter
  int outidx = 0;
  float filt_sum = 0;
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

// {{{ make_chirp_f
template <class T>
void make_chirp(T* full_in, int Nfull, T Fs, T Tfull) {
  T Ts = 1 / Fs;

  // Make input chirp
  for (int m = 0; m < Nfull; m++) {
    full_in[m] = sin(2 * M_PI * (Fs / Tfull * pow(Ts * m, 2) / 2));
  }
}
// }}}

