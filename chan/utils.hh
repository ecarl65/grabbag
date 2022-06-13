#include <iostream>
#include <complex>
#include <string>
#include <vector>
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
inline void write_out(std::string file, void *addr, size_t size, size_t numels) {
  FILE* fout = fopen(file.c_str(), "wb");
  fwrite(addr, size, numels, fout);
  fclose(fout);
}
// }}}

// {{{ write_append
inline void write_append(std::string file, void *addr, size_t size, size_t numels) {
  FILE* fout = fopen(file.c_str(), "ab");
  fwrite(addr, size, numels, fout);
  fclose(fout);
}
// }}}

// {{{ make_chirp
template <class T>
inline void make_chirp(std::vector<T> &full_in, T Fs, T Tfull) {
  T Ts = 1 / Fs;

  // Make input chirp
  for (size_t m = 0; m < full_in.size(); m++) {
    full_in[m] = sin(2 * M_PI * (Fs / Tfull * pow(Ts * m, 2) / 2));
  }
}
// }}}

