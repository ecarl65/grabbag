#!/usr/bin/env python3
import numpy as np
import pyudft
from scipy import signal
import matplotlib.pyplot as plt

# {{{ run_channelizer
def run_channelizer(downsample, sample_rate, filt_ord, full_ord, write=False, verbose=True):

    # Variables
    Nfilt = downsample * filt_ord + 1
    Nfull = 2**full_ord

    # Set up channelizer
    pu = pyudft.pudft(downsample, Nfilt, sample_rate, write, verbose)

    # Set up input data
    t = np.arange(Nfull) / sample_rate
    td = np.arange(len(t) / downsample) / (sample_rate / downsample)
    b = signal.chirp(t, 0, t[-1], 2 * sample_rate).astype('f')

    # Run channelizer
    od = pu.pyrun(b)

    # Display outputs
    print(f"Channelizer output dimensions: {od.shape}")
    plt.plot(t, b, color="dimgray", alpha=0.3)
    plt.plot(td, np.real(od[:, 0]))
    plt.plot(td, 2 * np.abs(od[:, 1:-1]))
    plt.plot(td, np.real(od[:, -1]))
    plt.title("Pybind11 UDFT")
    plt.xlabel("Time (sec)")
    plt.ylabel("Amplitude")
    plt.grid(True)

    plt.show()

# }}}

# {{{ __main__
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Test pybind11 Channelizer",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--downsample", help="Downsample amount", type=int, action="store",
            default=8)
    parser.add_argument("--sample-rate", help="Sample rate (Hz)", action="store",
            type=float, dest="sample_rate", default=10e3)
    parser.add_argument("--filt-ord", help="Filter order (length is this times downsample + 1)",
            action="store", type=int, default=8)
    parser.add_argument("--full-ord", help="Power of two for full number of samples",
            action="store", type=int, default=11)
    parser.add_argument("--write", help="Set to write output binary files",
            action="store_true")
    parser.add_argument("--verbose", help="Set to print verbose debugging",
            action="store_true")
    opts = parser.parse_args()

    run_channelizer(opts.downsample, opts.sample_rate, opts.filt_ord, opts.full_ord,
            opts.write, opts.verbose)

# }}}
