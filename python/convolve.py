#!/usr/bin/env python3

import numpy as np
import scipy as sp
import scipy.signal as sig
import matplotlib.pyplot as plt
import time

sample_rate = 100e3
N = 10000000
t = np.arange(N) / sample_rate
#  s = sig.chirp(t, 0, t[-1], sample_rate/2) * sig.gaussian(N)
#  s = sig.gausspulse(t - t[t.size >>1], fc=2e3, bw=0.001)
s = sig.gausspulse(t, fc=2e3, bw=0.001)

fig, axs = plt.subplots(1, 3)
axs[0].plot(t, s)

h = sig.firwin(128, 10e3, fs=sample_rate)

axs[1].plot(h)

print(sig.choose_conv_method(s, h))

tic1 = time.perf_counter()
s1 = sig.lfilter(h, 1, s)
toc1 = time.perf_counter()
tic2 = time.perf_counter()
s2 = sig.convolve(h, s)
#  s2 = sig.oaconvolve(h, s)
#  s2 = sig.fftconvolve(h, s)
toc2 = time.perf_counter()

axs[2].plot(s1[:10000])
axs[2].plot(s2[:10000], 'r.')

print(f"Length of s1: {len(s1)}")
print(f"Length of s2: {len(s2)}")
print(f"Time of lfilter: {toc1 - tic1}")
print(f"Time of convolve: {toc2 - tic2}")

plt.show()
