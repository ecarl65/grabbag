#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import signal
from pprint import pprint as pp
from shutil import get_terminal_size
np.set_printoptions(linewidth=get_terminal_size()[0])

# {{{ do_channelizer
def do_channelizer(downsample, oversample, sample_rate, filt, input_signal):
    """Do a naive version of the channelizer to test inputs and outputs against"""

    n_channels = downsample * oversample

    # Perform reshaping
    in_comm = np.reshape(
            np.r_[input_signal, np.zeros(n_channels - (len(input_signal) % n_channels))], 
            (n_channels, -1),
            order="F"
    )

    # Perform upsample
    up_comm = np.zeros((n_channels, in_comm.shape[1] * oversample))
    for chan in range(in_comm.shape[0]):
        up_comm[chan, ::oversample] = in_comm[chan, :]

    # Create filter
    filt_len = int(np.ceil(len(filt) / downsample) * downsample + n_channels)
    filt_cols = int(filt_len / downsample)
    filt_comm = np.zeros((n_channels, filt_cols))
    for rho in range(filt_comm.shape[0]):
        for ds_samp in range(filt_comm.shape[1]):
            idx = ds_samp * downsample - rho
            if 0 <= idx < len(filt):
                filt_comm[rho, ds_samp] = filt[idx]

    print("Filter:")
    pp(filt_comm)

    print("Data:")
    pp(up_comm[:, :10])

    # Do the filtering
    poly = np.zeros_like(up_comm)
    for ch, data_ch in enumerate(up_comm):
        #  poly[ch, :] = signal.lfilter(filt_comm[ch], 1, data_ch)
        poly[ch, :] = np.convolve(filt_comm[ch], data_ch, mode="same")

    # Do the vertical DFT
    channelized = np.fft.rfft(poly, axis=0)

    # plot_channelizer
    plot_channelizer(downsample, oversample, sample_rate, filt, input_signal, channelized.T)

# }}} 

# {{{ plot_channelizer
def plot_channelizer(downsample, oversample, sample_rate, filt, input_signal, channelized):

    n_channels = downsample * oversample

    rcParams["axes.grid"] = False
    fig, axs = plt.subplots(3, 2)
    axs[0, 0].plot(filt)
    axs[0, 0].set_title("Filter Taps")
    axs[0, 0].grid(True)

    fc_list = np.linspace(0, sample_rate/2, n_channels//2 + 1, endpoint=True)
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
    axs[0, 1].set_title("Input & Output (Normalized) \n (First/Last Real, Others Magnitude)")
    axs[0, 1].set_xlabel("Time (s)")
    axs[0, 1].set_ylabel("Amplitude")
    axs[0, 1].grid(True)

    mid_chan = (n_channels // 2 + 1) // 2
    tmp_data = channelized[:, mid_chan]
    axs[2, 0].plot(tds, np.real(tmp_data), label="real")
    axs[2, 0].plot(tds, np.imag(tmp_data), label="imag")
    axs[2, 0].plot(tds, np.abs(tmp_data), "g", label="mag")
    axs[2, 0].plot(tds, -np.abs(tmp_data), "g")
    axs[2, 0].grid(True)
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

    f, t, p = signal.stft(input_signal, fs=sample_rate, nperseg=n_channels)
    axs[1, 1].grid(False)
    axs[1, 1].pcolormesh(t, f, np.abs(p), shading="auto")
    axs[1, 1].set_title("STFT from SciPy")
    axs[1, 1].set_xlabel("Time")
    axs[1, 1].set_ylabel("Frequency")

    plt.tight_layout()

    plt.show()

# }}}

# {{{ main
def main(downsample, oversample, sample_rate, run):
    M = downsample
    I = oversample
    K = M * I
    Kh = K // 2 + 1

    indata = np.fromfile("input.bin", dtype=np.float32)
    outdata = np.fromfile("filtered.bin", dtype=np.float32)
    channelized = np.reshape(np.fromfile("channelized.bin", dtype=np.complex64), (-1, Kh))
    chann_buf = np.reshape(np.fromfile("onebuffer.bin", dtype=np.complex64), (-1, Kh))
    filt = np.fromfile("filter.bin", dtype=np.float32)
    fdata = np.reshape(np.fromfile("fftdata.bin", dtype=np.complex64), (K, -1))
    ffilt = np.reshape(np.fromfile("fftfilt.bin", dtype=np.complex64), (K, -1))

    if run:
        do_channelizer(downsample, oversample, sample_rate, filt, indata)
    else:
        plot_channelizer(downsample, oversample, sample_rate, filt, indata, channelized)

# }}}

# {{{ __main__
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plot channelizer output",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-d", "--downsample", help="Downsample amount", type=int, action="store",
            default=8)
    parser.add_argument("-o", "--oversample", help="Oversample factor", type=int, action="store",
            default=1)
    parser.add_argument("-s", "--sample-rate", help="Sample rate (Hz)", action="store",
            type=float, dest="sample_rate", default=10e3)
    parser.add_argument("--run", help="Run channelizer", action="store_true")
    opts = parser.parse_args()

    main(opts.downsample, opts.oversample, opts.sample_rate, opts.run)

# }}}
