#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

indata = np.reshape(np.fromfile("twodin.bin", dtype=np.float32), (32, -1));
outdata = np.reshape(np.fromfile("twodout.bin", dtype=np.complex64), (32, -1));

fs = 500
fig, axs = plt.subplots(2)
im0 = axs[0].pcolormesh(indata)
axs[0].set_title("Input Data")
f = (np.arange(outdata.shape[1] + 1) - 0.5) / outdata.shape[1] * fs
c = (np.arange(outdata.shape[0] + 1) - 0.5)
im1 = axs[1].pcolormesh(f, c, 20*np.log10(np.abs(outdata) + 1e-6))

axs[1].set_title("Rows FFT Output")
axs[1].set_xlabel("Frequency (index)")
axs[1].set_ylabel("Magnitude (dB)")

fig.colorbar(im0, ax=axs[0])
fig.colorbar(im1, ax=axs[1])

plt.tight_layout()

plt.show()
