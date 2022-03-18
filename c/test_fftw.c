/*
 * =====================================================================================
 *
 *       Filename:  test_fftw.c
 *
 *    Description: A test programm for FFTW 
 *
 *        Version:  1.0
 *        Created:  03/13/2022 06:35:07 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <tgmath.h>
#include <fftw3.h>
#include <stdint.h>

int main(int argc, char * argv[]) {
  // Declare variables
  const int N = 1024;
  const float pi = acosf(-1);
  float * vec = (float *) fftwf_malloc(sizeof(float) * N);
  fftwf_complex * output = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * N);

  // FFTW Plan
  int wise = fftwf_import_wisdom_from_filename("/home/eric/programming/c/rof1024");
  fftwf_plan p = fftwf_plan_dft_r2c_1d(N, vec, output, FFTW_EXHAUSTIVE);
  /* int writ = fftwf_export_wisdom_to_filename("/home/eric/programming/c/rof1024"); */

  // Populate input
  for (size_t idx = 0; idx < N; idx++) {
    vec[idx] = cosf(2 * pi * 0.05 * idx);
    if (idx < 10) {
      /* printf("%f\n", vec[idx]); */
    }
  }

  // Execute plan
  fftwf_execute(p); 

  // Output result
  for (size_t idx = 0; idx <= N/2; idx++) {
    /* printf("(%f+%fj)\n", creal(output[idx]), cimag(output[idx])); */
  }
  
  // Free memory
  fftwf_free(output);
  fftwf_free(vec);

  // 2 Dimensional FFT down one dimension only
  const int dim1 = 32;
  const int dim2 = 1024;
  float * vec2d = (float *) fftwf_malloc(sizeof(float) * dim1 * dim2);
  for (size_t idx1 = 0; idx1 < dim1; idx1++) {
    float freq = idx1 * 0.01f / 5.0;
    for (size_t idx2 = 0; idx2 < dim2; idx2++) {
      *(vec2d + idx1 * dim2 + idx2) = cosf(2 * pi * freq * idx2);
    }
  }
  fftwf_complex * ovec2d = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * dim1 * dim2);

  int n_size[] = {dim2};
  int rank = 1;
  int howmany = dim1;
  int idist = dim2;
  int odist = dim2;
  int istride = 1;
  int ostride = 1;
  int *inembed = n_size;
  int *onembed = n_size;
  fftwf_plan p2d = fftwf_plan_many_dft_r2c(rank, n_size, dim1, vec2d, inembed,
      istride, idist, ovec2d, onembed, ostride, odist, FFTW_ESTIMATE);
  /* int writ = fftwf_export_wisdom_to_filename("/home/eric/programming/c/rof1024"); */

  fftwf_execute(p2d);

  // Write output
  for (size_t idx1 = 0; idx1 < dim1; idx1++) {
    for (size_t idx2 = 0; idx2 < dim2; idx2++) {
      /* printf("%f ", *(vec2d + idx1 * dim2 + idx2)); */
      printf("(%f+%fj) ", creal(*(ovec2d + idx1 * dim2 + idx2)), cimag(*(ovec2d + idx1 * dim2 + idx2)));
    }
    printf("\n");
  }

  fftwf_free(vec2d);
  fftwf_free(ovec2d);

  return 0;
}

