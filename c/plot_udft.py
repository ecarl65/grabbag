#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import signal

rcParams["axes.grid"] = False

M = 8
Mh = int(M / 2) + 1

indata = np.fromfile("input.bin", dtype=np.float64)
outdata = np.fromfile("filtered.bin", dtype=np.float64)
channelized = np.reshape(np.fromfile("channelized.bin", dtype=np.complex128), (Mh, -1))
filt = np.fromfile("filter.bin", dtype=np.float64)
fdata = np.reshape(np.fromfile("fftdata.bin", dtype=np.complex128), (M, -1))
ffilt = np.reshape(np.fromfile("fftfilt.bin", dtype=np.complex128), (M, -1))

Nfilt = len(filt)
Noutdelay = int((Nfilt - 1) / (2 * M))

fs = 10e3
fig, axs = plt.subplots(3, 2)
axs[0, 0].plot(filt)
axs[0, 0].set_title("Filter Taps")
axs[0, 0].grid(True)

tin = np.arange(len(indata)) / fs
tout = np.arange(len(outdata)) / fs
w, h = signal.freqz(filt, 1, fs=10e3)
axs[1, 0].plot(w, 20*np.log10(np.abs(h)))
axs[1, 0].plot([1250/2, 1250/2], [-100, 0], "r:")
axs[1, 0].set_title("Filter Frequency Response")
axs[1, 0].set_xlabel("Frequency (Hz)")
axs[1, 0].set_ylabel("Magnitude (dB)")
axs[1, 0].grid(True)

pchan = 0
tds = np.arange(channelized.shape[1]) / fs * M
axs[0, 1].plot(tin, indata, label="Input", color="dimgray", alpha=0.3)
axs[0, 1].plot(tds[:(len(tds) - Noutdelay)], np.real(channelized[0, Noutdelay:].T) / channelized.shape[1], label=f"Channel 0 (real)")
axs[0, 1].plot(tds[:(len(tds) - Noutdelay)], np.abs(channelized[1:-1, Noutdelay:].T) / channelized.shape[1] * 2, label="Other Channels (mag)")
axs[0, 1].plot(tds[:(len(tds) - Noutdelay)], np.real(channelized[-1, Noutdelay:].T) / channelized.shape[1], label=f"Channel {M >> 1} (real)")
axs[0, 1].set_title("Input & Output (Normalized)")
axs[0, 1].set_xlabel("Time (s)")
axs[0, 1].set_ylabel("Amplitude")
axs[0, 1].legend()
axs[0, 1].grid(True)

f = np.linspace(0, fs / 2, fdata.shape[1], endpoint=False)
c = np.arange(fdata.shape[0])
axs[2, 0].pcolormesh(f, c, np.abs(fdata), shading="auto")
axs[2, 0].set_title("FFT of Polyphase Data")
axs[2, 0].set_xlabel("Frequency (Hz)")
axs[2, 0].set_ylabel("Polyphase Channel")

ch_nd = channelized[:, Noutdelay:]
c = np.arange(ch_nd.shape[0])
t = np.arange(ch_nd.shape[1]) / fs * M
axs[2, 1].pcolormesh(t, c, np.abs(ch_nd), shading="auto")
axs[2, 1].set_title("Channelized Output")
axs[2, 1].set_ylabel("Channel")
axs[2, 1].set_xlabel("Time Sample")

f, t, p = signal.stft(indata, fs=10e3, nperseg=M)
axs[1, 1].pcolor(t, f, np.abs(p), shading="auto")
axs[1, 1].set_title("STFT from SciPy")
axs[1, 1].set_xlabel("Time")
axs[1, 1].set_ylabel("Frequency")

plt.tight_layout()

plt.show()
