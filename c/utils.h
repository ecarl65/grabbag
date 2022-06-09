#include <stdlib.h>
#include <stdint.h>
#include <tgmath.h>
#include <string.h>
#include <time.h>

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
void write_out(char *file, void *addr, size_t size, size_t numels) {
  FILE* fout = fopen(file, "ab");
  fwrite(addr, size, numels, fout);
  fclose(fout);
}
// }}}

// {{{ coarse_gaussian_seed
double coarse_gaussian_seed(unsigned int seed) {
  srand(seed);
  double out = 0;
  for (int m = 0; m < 12; m++) {
    out += (double) (rand() % 100);
  }
  out /= 100.0;
  out -= 6.0;

  return out;
}
// }}}

// {{{ coarse_gaussian
double coarse_gaussian() {
  unsigned int seed = time(NULL);
  double output = coarse_gaussian_seed(seed);

  return output;
}
// }}}

// {{{ poly_filt_design
void poly_filt_design(int Nfilt, double fcutoff, double Fs, double* filt, double* filt_full, int ncols, int downsamp) {
  double Ts = 1 / Fs;
  const int Nfilt_half = (Nfilt - 1) >> 1;

  // Make filter
  int outidx = 0;
  double filt_sum = 0;
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
}
// }}} 

// {{{ make_chirp
void make_chirp(double* full_in, int Nfull, double Fs, double Tfull) {
  double Ts = 1 / Fs;

  // Make input chirp
  for (int m = 0; m < Nfull; m++) {
    // full_in[m] = sin(2 * M_PI * (Fs / Tfull * pow(Ts * m, 2) / 2)) + coarse_gaussian() / 10.0;
    full_in[m] = sin(2 * M_PI * (Fs / Tfull * pow(Ts * m, 2) / 2));
  }
}
// }}}

// {{{ make_sinusoid
void make_sinuosoid(double* full_in, int Nfull, double Fs, double fc) {
  double Ts = 1 / Fs;

  // Make input sinusoid
  for (int m = 0; m < Nfull; m++) {
    full_in[m] = cos(2 * M_PI * fc * m * Ts) + coarse_gaussian() / 10.0;
  }
}
// }}}
