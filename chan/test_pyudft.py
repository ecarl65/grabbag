#!/usr/bin/env python3
"""Run the pybind11 channelizer"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import time
import pyudft
from plot_udft import plot_channelizer

# {{{ run_channelizer
def run_channelizer(downsample, oversample, sample_rate, filt_ord, full_ord, write=False, verbose=True, no_plots=False):
    """Main entry point to call the channelizer"""

    # Variables
    n_channels = downsample * oversample
    n_filt = downsample * oversample * filt_ord + 1
    n_full = 2**full_ord

    # Set up channelizer
    pu = pyudft.pudft(downsample, oversample, n_filt, sample_rate, write, verbose)

    # Get filter
    filt = pu.get_filter()

    # Set up input data
    time_vec = np.arange(n_full) / sample_rate
    time_vec_ds = np.arange(len(time_vec) / downsample) / (sample_rate / downsample)
    input_signal = signal.chirp(time_vec, 0, time_vec[-1], 2 * sample_rate).astype('f')

    # Run channelizer
    start = time.time()
    channelized = pu.pyrun(input_signal)
    stop = time.time()
    print(f"Elapsed time to run: {stop - start}")

    # Display outputs
    if verbose:
        print(f"Channelizer output dimensions: {channelized.shape}")

    # Plot dependencies: filt, sample_rate, downsample, input_signal, channelized
    if not no_plots:
        plot_channelizer(downsample, oversample, sample_rate, filt, input_signal, channelized)

# }}}

# {{{ __main__
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Test pybind11 Channelizer",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-d", "--downsample", help="Downsample amount", type=int, action="store",
            default=8)
    parser.add_argument("-o", "--oversample", help="Oversample factor", type=int, action="store",
            default=1)
    parser.add_argument("-s", "--sample-rate", help="Sample rate (Hz)", action="store",
            type=float, dest="sample_rate", default=10e3)
    parser.add_argument("--filt-ord", help="Filter order (length is this times downsample + 1)",
            action="store", type=int, default=8)
    parser.add_argument("--full-ord", help="Power of two for full number of samples",
            action="store", type=int, default=11)
    parser.add_argument("-w", "--write", help="Set to write output binary files",
            action="store_true")
    parser.add_argument("-v", "--verbose", help="Set to print verbose debugging",
            action="store_true")
    parser.add_argument("-n", "--no-plots", help="Set to suppress plots", action="store_true")
    opts = parser.parse_args()

    run_channelizer(opts.downsample, opts.oversample, opts.sample_rate, opts.filt_ord, opts.full_ord,
            opts.write, opts.verbose, opts.no_plots)

# }}}
