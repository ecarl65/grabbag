#!/usr/bin/env python

import numpy as np
import scipy as sp
from scipy import signal
import matplotlib.pyplot as plt
from copy import deepcopy

fs = 25
npoints = int(500)

fc = 0.2
t = np.arange(npoints) / fs
s = np.cos(2 * np.pi * fc * t)

h_len = 11
#  h = deepcopy(s[:128])
#  h /= 64
h = np.ones(h_len) / h_len

w = np.random.randn(len(t)) / 20

s_n = s + w

# numpy
cf = np.convolve(s_n, h, mode='full')
cv = np.convolve(s_n, h, mode='valid')
cs = np.convolve(s_n, h, mode='same')

# fft approach
#  full_len = len(h) + len(s_n) - 1
full_len = 2 * npoints
ff = np.fft.ifft(np.conj(np.fft.fft(h, full_len)) * np.fft.fft(s_n, full_len))
t_ff = np.arange(len(ff)) / fs

# overlap save
nos = 4 * h_len
nv = nos - h_len + 1  # Not sure about plus 1, must verify
n_buffs = int(npoints / nv)
H = np.fft.fft(h, nos)
sfo = np.zeros(len(s_n))
for x in range(n_buffs):
    st = nv * x
    sp = st + nos
    D = np.fft.fft(s_n[st:sp], nos)
    do = np.fft.ifft(H * D)
    dov = do[nos-nv:]
    sfo[st:st + nv] = np.real(dov)

tsfo = np.arange(len(sfo)) / fs

delay_v = ((h_len - 1) // 2) / fs

# Test
#  np.allclose(sfo[:959], cv[:959])

tf = np.arange(len(cf)) / fs
plt.plot(t, s, label='orig')
plt.plot(t, s_n, label='noisy')
#  plt.plot(t[:len(h)], h, label='filter')
#  plt.plot(tf - delay_v, cf, label='full', linestyle='-.')
#  plt.plot(tf[:len(cv)] + delay_v, cv, label='valid', linestyle=':')
#  plt.plot(tf[:len(cs)], cs, label='same')
#  plt.plot(t_ff - delay_v, np.real(ff), label='ov fft', color='magenta', marker='o', linestyle='None')
plt.plot(t_ff, (np.real(ff)), label='ov fft', color='magenta')
#  plt.plot(tsfo + delay_v, np.real(sfo), label='o/s', marker='+', linestyle='None')
plt.legend()
plt.show()
