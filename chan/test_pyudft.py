#!/usr/bin/env python3
### current bpython session - make changes and save to reevaluate session.
### lines beginning with ### will be ignored.
### To return to bpython without reevaluating make no changes to this file
### or save an empty file.
import numpy as np
import pyudft
from scipy import signal
import matplotlib.pyplot as plt

M = 8
Nfilt = M * 8 + 1
Fs = 10e3
write_out = True
verbose = True

pu = pyudft.pudft(M, Nfilt, Fs, write_out, verbose)

Nfull = 2048
t = np.arange(Nfull) / Fs
b = signal.chirp(t, 0, t[-1], 2 * Fs).astype('f')
od = pu.pyrun(b)
td = np.arange(len(t) / M) / (Fs / M)

plt.plot(t, b, color="dimgray", alpha=0.2)
plt.plot(td, np.real(od[:, 0]))
plt.plot(td, 2 * np.abs(od[:, 1:-1]))
plt.plot(td, np.real(od[:, -1]))
plt.title("Pybind11 UDFT")
plt.xlabel("Time (sec)")
plt.ylabel("Amplitude")

plt.show()
