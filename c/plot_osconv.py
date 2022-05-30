#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

indata = np.fromfile("input.bin", dtype=np.float64);
outdata = np.fromfile("filtered.bin", dtype=np.float64);
filt = np.fromfile("filter.bin", dtype=np.float64);

fs = 10e3
fig, axs = plt.subplots(3)
axs[0].plot(filt)
axs[0].set_title("Filter Taps")
tin = np.arange(len(indata)) / fs
tout = np.arange(len(outdata)) / fs
w, h = signal.freqz(filt, 1, fs=10e3)
axs[1].plot(w, 20*np.log10(np.abs(h)))
axs[1].plot([100, 100], [-60, 0], "r:")
axs[1].set_title("Filter Frequency Response")
axs[1].set_xlabel("Frequency (Hz)")
axs[1].set_ylabel("Magnitude (dB)")
axs[2].plot(tin, indata, label="Input")
axs[2].plot(tout, outdata, label="O/S Filtered");
axs[2].set_title("Input & Output")
axs[2].set_xlabel("Time (s)")
axs[2].set_ylabel("Amplitude")
axs[2].legend()

plt.show()
