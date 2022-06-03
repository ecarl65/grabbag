#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

indata = np.reshape(np.fromfile("twodin.bin", dtype=np.float32), (32, -1))
outdata = np.reshape(np.fromfile("twodout.bin", dtype=np.complex64), (32, -1))
oned = np.fromfile("oned.bin", dtype=np.complex64)

plt.rcParams["axes.grid"] = False

fs = 500
fig, axs = plt.subplots(2, 2)
f1 = np.linspace(0, 0.5, len(oned), endpoint=False)
axs[0, 0].plot(f1, 20 * np.log10(np.abs(oned)))
axs[0, 0].set_title("0.05 Peak")
axs[0, 0].set_xlabel("Norm. Freq")
axs[0, 0].set_ylabel("Magnitude")
axs[0, 0].grid(True)

im0 = axs[0, 1].pcolormesh(indata)
axs[0, 1].set_title("Input Data")

f = np.linspace(0, fs / 2, outdata.shape[1], endpoint=False)
c = np.arange(outdata.shape[0])
im1 = axs[1, 1].pcolormesh(f, c, 20*np.log10(np.abs(outdata) + 1e-6), shading="nearest")

axs[1, 1].set_title("Rows FFT Output")
axs[1, 1].set_xlabel("Frequency (Hz)")
axs[1, 1].set_ylabel("Magnitude (dB)")

axs[1, 0].plot(f, np.abs(outdata[-1, :]))
axs[1, 0].set_xlabel("Frequency (Hz)")
axs[1, 0].set_ylabel("FFT Output Magnitude")
axs[1, 0].set_title("15.5 Hz Row")
axs[1, 0].grid(True)

fig.colorbar(im0, ax=axs[0, 1])
fig.colorbar(im1, ax=axs[1, 1])

plt.tight_layout()
plt.rcParams["axes.grid"] = True

plt.show()
