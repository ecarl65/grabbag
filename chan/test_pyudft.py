#!/usr/bin/env python3
"""Run the pybind11 channelizer"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import pyudft

# {{{ run_channelizer
def run_channelizer(downsample, sample_rate, filt_ord, full_ord, write=False, verbose=True):
    """Main entry point to call the channelizer"""

    # Variables
    n_filt = downsample * filt_ord + 1
    n_full = 2**full_ord

    # Set up channelizer
    pu = pyudft.pudft(downsample, n_filt, sample_rate, write, verbose)

    # Get filter
    filt = pu.get_filter()

    # Set up input data
    time_vec = np.arange(n_full) / sample_rate
    time_vec_ds = np.arange(len(time_vec) / downsample) / (sample_rate / downsample)
    input_signal = signal.chirp(time_vec, 0, time_vec[-1], 2 * sample_rate).astype('f')

    # Run channelizer
    channelized = pu.pyrun(input_signal)

    # Display outputs
    if verbose:
        print(f"Channelizer output dimensions: {channelized.shape}")

    fig, axs = plt.subplots(3, 2)
    axs[0, 0].plot(filt)
    axs[0, 0].set_title("Filter Taps")
    axs[0, 0].grid(True)

    fc_list = np.linspace(0, sample_rate/2, downsample//2 + 1, endpoint=True)
    tin = np.arange(len(input_signal)) / sample_rate
    tfilt = np.arange(len(filt)) / sample_rate

    for fc in fc_list:
        w, h = signal.freqz(filt * np.exp(1j * 2 * np.pi * fc * tfilt), 1, fs=sample_rate, whole=False)
        axs[1, 0].plot(w, 20*np.log10(np.abs(h)))
        axs[1, 0].plot([1250/2, 1250/2], [-100, 0], "r:")
        axs[1, 0].plot([fc, fc], [-100, 0], "m:")
    axs[1, 0].set_title("Channelizer Frequency Response")
    axs[1, 0].set_xlabel("Frequency (Hz)")
    axs[1, 0].set_ylabel("Magnitude (dB)")
    axs[1, 0].grid(True)

    tds = np.arange(channelized.shape[0]) / sample_rate * downsample
    axs[0, 1].plot(tin, input_signal, label="Input", color="dimgray", alpha=0.3)
    axs[0, 1].plot(tds, np.real(channelized[:, 0]), label=f"Channel 0 (real)")
    axs[0, 1].plot(tds, np.abs(channelized[:, 1:-1]) * 2, label="Other Channels (mag)")
    axs[0, 1].plot(tds, np.real(channelized[:, -1]), label=f"Channel {downsample >> 1} (real)")
    axs[0, 1].set_title("Input & Output (Normalized)")
    axs[0, 1].set_xlabel("Time (s)")
    axs[0, 1].set_ylabel("Amplitude")
    axs[0, 1].text(0.1, -1.5, "First/Last Channels real, others magnitude", horizontalalignment="center", verticalalignment="center", bbox=dict(facecolor='white', alpha=0.9))
    axs[0, 1].grid(True)
    axs[0, 1].set_ylim(-2, 1.5)

    mid_chan = (downsample // 2 + 1) // 2
    tmp_data = channelized[:, mid_chan]
    axs[2, 0].plot(tds, np.real(tmp_data), label="real")
    axs[2, 0].plot(tds, np.imag(tmp_data), label="imag")
    axs[2, 0].plot(tds, np.abs(tmp_data), "g", label="mag")
    axs[2, 0].plot(tds, -np.abs(tmp_data), "g")
    axs[2, 0].set_title(f"Output Channel {mid_chan}")
    axs[2, 0].set_xlabel("Time (sec)")
    axs[2, 0].set_ylabel("Amplitude")
    axs[2, 0].legend()

    c = np.arange(channelized.shape[1])
    t = np.arange(channelized.shape[0]) / sample_rate * downsample
    axs[2, 1].grid(False)
    axs[2, 1].pcolormesh(t, c, np.abs(channelized.T), shading="auto")
    axs[2, 1].set_title("Channelized Output")
    axs[2, 1].set_ylabel("Channel")
    axs[2, 1].set_xlabel("Time Sample")

    f, t, p = signal.stft(input_signal, fs=sample_rate, nperseg=downsample)
    axs[1, 1].grid(False)
    axs[1, 1].pcolormesh(t, f, np.abs(p), shading="auto")
    axs[1, 1].set_title("STFT from SciPy")
    axs[1, 1].set_xlabel("Time")
    axs[1, 1].set_ylabel("Frequency")

    plt.tight_layout()

    plt.show()
# }}}

# {{{ __main__
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Test pybind11 Channelizer",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-d", "--downsample", help="Downsample amount", type=int, action="store",
            default=8)
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
    opts = parser.parse_args()

    run_channelizer(opts.downsample, opts.sample_rate, opts.filt_ord, opts.full_ord,
            opts.write, opts.verbose)

# }}}
