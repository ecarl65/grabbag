#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

M = 8
Mh = int(M / 2) + 1

indata = np.fromfile("input.bin", dtype=np.float64)
outdata = np.fromfile("filtered.bin", dtype=np.float64)
channelized = np.reshape(np.fromfile("channelized.bin", dtype=np.complex128), (Mh, -1))
filt = np.fromfile("filter.bin", dtype=np.float64)
fdata = np.reshape(np.fromfile("fftdata.bin", dtype=np.complex128), (M, -1))
ffilt = np.reshape(np.fromfile("fftfilt.bin", dtype=np.complex128), (M, -1))

fs = 10e3
fig, axs = plt.subplots(3, 2)
axs[0, 0].plot(filt)
axs[0, 0].set_title("Filter Taps")
tin = np.arange(len(indata)) / fs
tout = np.arange(len(outdata)) / fs
w, h = signal.freqz(filt, 1, fs=10e3)
axs[1, 0].plot(w, 20*np.log10(np.abs(h)))
axs[1, 0].plot([100, 100], [-60, 0], "r:")
axs[1, 0].set_title("Filter Frequency Response")
axs[1, 0].set_xlabel("Frequency (Hz)")
axs[1, 0].set_ylabel("Magnitude (dB)")

pchan = 1
tds = np.arange(channelized.shape[1]) / fs * M
axs[0, 1].plot(tin, indata, label="Input", color="dimgray", alpha=0.3)
axs[0, 1].plot(tds, np.real(channelized[pchan, :]) / channelized.shape[1], label=f"Channel {pchan} (real)")
axs[0, 1].plot(tds, np.imag(channelized[pchan, :]) / channelized.shape[1], label=f"Channel {pchan} (imag)")
axs[0, 1].set_title("Input & Output")
axs[0, 1].set_xlabel("Time (s)")
axs[0, 1].set_ylabel("Amplitude")
axs[0, 1].legend()

#  axs[1, 1].pcolormesh(np.abs(ffilt))
#  axs[1, 1].set_title("Filter Decomposition")
#  axs[1, 1].set_xlabel("Frequency Index")
#  axs[1, 1].set_ylabel("Polyphase Channel")

axs[2, 1].pcolormesh(np.abs(fdata))
axs[2, 1].set_title("FFT of Polyphase Data")
axs[2, 1].set_xlabel("Frequency Index")
axs[2, 1].set_ylabel("Polyphase Channel")

axs[2, 0].pcolormesh(np.abs(channelized))
axs[2, 0].set_title("Channelized Output")
axs[2, 0].set_ylabel("Channel")
axs[2, 0].set_xlabel("Time Sample")

f, t, p = signal.stft(indata, fs=10e3, nperseg=M)
axs[1, 1].pcolormesh(t, f, np.abs(p))
axs[1, 1].set_title("STFT from SciPy")
axs[1, 1].set_xlabel("Time")
axs[1, 1].set_ylabel("Frequency")

plt.tight_layout()

plt.show()
