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

class UDFT {
  public:
    // Variables
    int downsamp, n_full, n_filt;
    float samp_rate;
    int n_buffer, n_cols, n_cols_fft, n_rows_fft, n_fft_h, n_fft_v, n_out, n_delay, n_delay_r, n_delay_samp;
    int idx_out_valid_r, idx_out_valid_samp, n_in_valid, n_out_valid;
    float samp_period, f_cutoff;
    struct fft_config fwd_c, filt_c, inv_c, col_c;
    float *filt, *buffer_in, *filt_full, *conv_out;
    std::complex<float> *full_out, *fft_in, *fft_filt, *fft_mult, *udft;
    fftwf_plan psig, pinv, pudft;
    bool write, debug;

    void poly_filt_design();

    UDFT(int downsamp, int n_full, int n_filt, float samp_rate, bool write = false, bool debug = true);

    ~UDFT();

    void run(float *indata);

};

