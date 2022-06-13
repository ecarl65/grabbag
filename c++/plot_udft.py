#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import signal

rcParams["axes.grid"] = False

if len(sys.argv) > 1:
    M = int(sys.argv[1])
else:
    M = 8
Mh = int(M / 2) + 1

indata = np.fromfile("input.bin", dtype=np.float32)
outdata = np.fromfile("filtered.bin", dtype=np.float32)
channelized = np.reshape(np.fromfile("channelized.bin", dtype=np.complex64), (-1, Mh))
chann_buf = np.reshape(np.fromfile("onebuffer.bin", dtype=np.complex64), (-1, Mh))
filt = np.fromfile("filter.bin", dtype=np.float32)
fdata = np.reshape(np.fromfile("fftdata.bin", dtype=np.complex64), (M, -1))
ffilt = np.reshape(np.fromfile("fftfilt.bin", dtype=np.complex64), (M, -1))

Nfilt = len(filt)
Noutdelay = int((Nfilt - 1) / (2 * M))

fs = 10e3
fig, axs = plt.subplots(3, 2)
axs[0, 0].plot(filt)
axs[0, 0].set_title("Filter Taps")
axs[0, 0].grid(True)

fc_list = np.linspace(0, fs/2, M//2 + 1, endpoint=True)
tin = np.arange(len(indata)) / fs
tout = np.arange(len(outdata)) / fs
tfilt = np.arange(len(filt)) / fs

for fc in fc_list:
    w, h = signal.freqz(filt * np.exp(1j * 2 * np.pi * fc * tfilt), 1, fs=10e3, whole=False)
    axs[1, 0].plot(w, 20*np.log10(np.abs(h)))
    axs[1, 0].plot([1250/2, 1250/2], [-100, 0], "r:")
    axs[1, 0].plot([fc, fc], [-100, 0], "m:")
axs[1, 0].set_title("Channelizer Frequency Response")
axs[1, 0].set_xlabel("Frequency (Hz)")
axs[1, 0].set_ylabel("Magnitude (dB)")
axs[1, 0].grid(True)

pchan = 0
tds = np.arange(channelized.shape[0]) / fs * M
tob = np.arange(chann_buf.shape[0]) / fs * M
axs[0, 1].plot(tin, indata, label="Input", color="dimgray", alpha=0.3)
axs[0, 1].plot(tds, np.real(channelized[:, 0]), label=f"Channel 0 (real)")
axs[0, 1].plot(tds, np.abs(channelized[:, 1:-1]) * 2, label="Other Channels (mag)")
#  axs[0, 1].plot(tds, np.abs(channelized[:, 1:-1]), label="Other Channels (mag)")
axs[0, 1].plot(tds, np.real(channelized[:, -1]), label=f"Channel {M >> 1} (real)")
#  axs[0, 1].plot(tob, np.abs(chann_buf) / chann_buf.shape[0] * 2, label=f"One Buffer Channelized Output")
axs[0, 1].set_title("Input & Output (Normalized)")
axs[0, 1].set_xlabel("Time (s)")
axs[0, 1].set_ylabel("Amplitude")
#  axs[0, 1].legend()
axs[0, 1].text(0.1, -1.5, "First/Last Channels real, others magnitude", horizontalalignment="center", verticalalignment="center", bbox=dict(facecolor='white', alpha=0.9))
axs[0, 1].grid(True)
axs[0, 1].set_ylim(-2, 1.5)

f = np.linspace(0, fs / 2, fdata.shape[1], endpoint=False)
c = np.arange(fdata.shape[0])
axs[2, 0].pcolormesh(f, c, np.abs(fdata), shading="auto")
axs[2, 0].set_title("FFT of Polyphase Data")
axs[2, 0].set_xlabel("Frequency (Hz)")
axs[2, 0].set_ylabel("Polyphase Channel")

c = np.arange(channelized.shape[1])
t = np.arange(channelized.shape[0]) / fs * M
axs[2, 1].pcolormesh(t, c, np.abs(channelized.T), shading="auto")
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
