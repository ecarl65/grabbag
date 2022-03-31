#!/usr/bin/env python3

import bluefile as bf
import threedb
import threedb.xmalt as xm
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

fs = 2.5e3
dur = 1.0
tdoa_range = 20e-3
fdoa_range = 100.0
points_per_sym = 20

num_points = int(dur * fs / points_per_sym) * points_per_sym
tight = 4
snr_min = 50

num_signals = 3

# Just one with ambiguities

sig1 = np.zeros(num_points, dtype=np.complex128)
sig2 = np.zeros(num_points, dtype=np.complex128)
for sig_num in range(num_signals):
    approx_snr = np.random.random() * 20 + snr_min
    amp = 10 ** (approx_snr / 20)

    tdoa = (np.random.random() * tdoa_range - tdoa_range / 2) / 2
    fdoa = (np.random.random() * fdoa_range - fdoa_range / 2) / 2

    tdoa_samps = int(tdoa * fs)

    data = amp * signal.lfilter(np.ones(10), 10, threedb.cnormal(scale=1, size=num_points))
    sig1 += data

    data_mod = np.roll(data, tdoa_samps)
    t = np.arange(num_points) / fs
    data_mod = data_mod * np.exp(1j * 2 * np.pi * fdoa * t)

    sig2 += data_mod

    noise1 = threedb.cnormal(scale=1, size=num_points)
    noise2 = threedb.cnormal(scale=1, size=num_points)

    sig1 += noise1
    sig2 += noise2

    print(f"amplitude: {amp}, tdoa: {tdoa * 1e3} ms, fdoa: {fdoa} Hz")

hdr = bf.header(type=1000, format="CD", xstart=0, xunits=1, xdelta=1/fs)
bf.write('signal1.tmp', hdr=hdr, data=sig1)
bf.write('signal2.tmp', hdr=hdr, data=sig2)

print(f"TDOA Range: {tdoa_range}")
print(f"FDOA Range: {fdoa_range}")
