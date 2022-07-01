#!/usr/bin/env python3
"""Test the Gain from FFT/PSD calculations. Why is there a discrepancy. For the math see:
https://ccrma.stanford.edu/~jos/sasp/Processing_Gain.html
"""

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.mlab import psd

def cnormal(N):
    return np.dot(np.random.randn(N, 2), [1 / np.sqrt(2), 1j / np.sqrt(2)])

def a2db(x):
    return 20 * np.log10(np.abs(x))

def p2db(x):
    return 10 * np.log10(np.abs(x))

# Constants
num_samps = 8192
samp_rate = 2.5e6
snr = 20
nfft = 1024
pad = 4096

# Input signal and noise
time = np.arange(num_samps) / samp_rate
freq_0 = samp_rate / nfft * 10.33  # Alter fractional part to hit bin center or not
complex_noise = cnormal(len(time)) * 10 ** (-snr / 20)
complex_signal = np.exp(1j * 2 * np.pi * freq_0 * time)
# complex_signal = 1

# To normalize all outputs the same treat FFT like multiplying by rect window and normalize "filter"
h_rect = np.ones(nfft) / nfft
spn = complex_signal + complex_noise
SPN = np.fft.fftshift(np.fft.fft(h_rect * spn[:nfft], pad))

# My own periodogram. Normalize filter as well for same peak power. No overlap.
win = signal.windows.hamming(nfft)
win /= np.sum(win)
n_loops = int(num_samps / nfft)
SPN_per = np.zeros(pad)
for iter in range(n_loops):
    SPN_per += np.abs(np.fft.fftshift(np.fft.fft(win * spn[iter*nfft:(iter+1)*nfft], pad)))**2
SPN_per /= n_loops

# Matlab PSD
Pxx, freqs = psd(spn, Fs=samp_rate, NFFT=nfft, sides="twosided", noverlap=0, pad_to=pad, scale_by_freq=False)

# Plot results
plt.close("all")
fig = plt.figure(figsize=(24,12), tight_layout=True)
axs = fig.add_subplot(111)
axs.plot(freqs, a2db(SPN), color="b", label="FFT", alpha=0.5)
axs.plot(freqs, p2db(SPN_per), color="r", label="Periodogram", alpha=0.5)
axs.plot(freqs, p2db(Pxx), color="k", label="PSD")
axs.set_title("PSD of Sinusoid in Noise")
axs.set_xlabel("Frequency (Hz)")
axs.set_ylabel("Normalized Power (dB)")
axs.legend()
plt.show()

# Better to do both sides but for now just get noise stats on one side of signal
ncutoff = int(pad * 0.4)

# Expected gains and SNRs
proc_gain = p2db(nfft)
expect = proc_gain + snr

# Print out some results
print(f"Input SNR: {snr}")

print(f"Processing gain: {proc_gain:.2f}")
print(f"Expected SNR: {expect:.2f}")
print(f"Number of averages in periodogram: {n_loops}")

# Quantify SNR per method
noise_mean = {
    "PSD": np.mean(Pxx[:ncutoff]),
    "FFT": np.mean(np.abs(SPN[:ncutoff])**2),
    "Per": np.mean(SPN_per[:ncutoff])
}
for psd_type, nm in noise_mean.items():
    snr_db = p2db((1 - nm) / nm)
    # snr_db = p2db(1 / nm)
    print(f"{psd_type}: SNR = {snr_db:.2f}")
    print(f"{psd_type}: Delta = {expect - snr_db:.2f}")

