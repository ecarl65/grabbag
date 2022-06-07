#include <stdlib.h>
#include <stdint.h>
#include <complex.h>
#include <tgmath.h>
#include <string.h>

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

void write_out(char *file, void *addr, size_t size, size_t numels) {
  FILE* fout = fopen(file, "wb");
  fwrite(addr, size, numels, fout);
  fclose(fout);
}

double coarse_gaussian() {
  double out = 0;
  for (int m = 0; m < 12; m++) {
    out += (double) (rand() % 100);
  }
  out /= 100.0;
  out -= 6.0;
}

