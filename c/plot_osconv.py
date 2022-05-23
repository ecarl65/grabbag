#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

indata = np.fromfile("input.bin", dtype=np.float64);
outdata = np.fromfile("filtered.bin", dtype=np.float64);
filt = np.fromfile("filter.bin", dtype=np.float64);

fig, axs = plt.subplots(2)
axs[0].plot(filt)
tin = np.arange(len(indata))
tout = np.arange(len(outdata))
axs[1].plot(tin, indata)
axs[1].plot(tout, outdata);

plt.show()
