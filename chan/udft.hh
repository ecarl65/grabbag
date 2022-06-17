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
#include "utils.hh"

using cfloat = std::complex<float>;

class UDFT {
  public:

    UDFT(int downsamp, int oversamp, int n_filt, float samp_rate, bool write = false, bool debug = true);
    ~UDFT();
    void poly_filt_design();
    std::vector<std::vector<cfloat>> run(float *indata, int n_full);

    // Variables
    int downsamp, oversamp, n_channels, n_filt;
    float samp_rate;
    int n_buffer, n_cols_filt, n_cols_data, n_cols_filt_fft, n_cols_data_fft, n_channels_out;
    int n_tot_fft_data, n_tot_fft_filt, n_tot_out;
    int n_delay, n_delay_r, n_delay_samp;
    int idx_out_valid_r, idx_out_valid_samp, n_in_valid, n_out_valid_r, n_out_valid_samp;
    float samp_period, f_cutoff;
    fft_config fwd_c, filt_c, inv_c, col_c;
    float *filt, *buffer_in, *filt_full, *conv_out;  // FFTW real arrays
    cfloat *fft_in, *fft_filt, *fft_mult, *udft;     // FFTW complex arrays
    fftwf_plan psig, pinv, pudft;
    bool write, debug;
};

