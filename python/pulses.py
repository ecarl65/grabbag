#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from pprint import pprint

matplotlib.use("TkAgg")


def main(opts):
    """Main Call"""

    padding = opts.duration
    total_dur = 2 * padding + opts.duration

    # TODO Try a variety of leading edge durations
    tr = opts.rise_percent / 100 * opts.duration

    print(f"Rising edge time: {tr:.3e}")

    t_0 = padding
    t_1 = padding + opts.duration
    le_st = t_0 - tr / 2
    le_sp = t_0 + tr / 2
    te_st = t_1 - tr / 2
    te_sp = t_1 + tr / 2

    t = np.arange(0, total_dur, 1 / opts.sample_rate)

    # Set up zeros and ones portion
    p = np.zeros(t.size)
    p[(t > le_sp) & (t < te_st)] = 1

    # Set up leading and trailing edges
    le_idx = (t >= le_st) & (t <= le_sp)
    te_idx = (t >= te_st) & (t <= te_sp)

    if opts.type == 'lin':
        # Leading trailing edges for linear
        b = -le_st / tr
        p[le_idx] = t[le_idx] / tr - le_st / tr
        p[te_idx] = -t[te_idx] / tr + te_sp / tr
    elif opts.type == 'sin':
        p[le_idx] = 0.5 * (1 + np.sin(np.pi / tr * (t[le_idx] - t_0)))
        p[te_idx] = 0.5 * (1 - np.sin(np.pi / tr * (t[te_idx] - t_1)))

    compute_diff(t, p, opts, tr)

    return t, p


def plot_freq(t, pulse, opts):
    """Plot the frequency and compute G4"""

    fft_size = 2 ** int(np.ceil(np.log2(pulse.size * 10)))

    Pxx = np.fft.fft(pulse, fft_size)

    df = opts.sample_rate / fft_size
    freqs = np.arange(0, opts.sample_rate, df)
    denom = np.sum(np.abs(Pxx) ** 2 * df)
    num = np.sum(np.abs(Pxx) ** 2 * df * freqs ** 2)

    G4 = num / denom
    print(f"Spectral mean squared spread: {G4:.3e}")

    fig, axs = plt.subplots(2, 1)
    axs[0].plot(t * 1e6, p)
    axs[0].set_title(f'Pulse (G4 = {G4:.3e})')
    axs[0].set_ylabel('Amplitude')
    axs[0].set_xlabel(r'Time $(\mu s)$')

    axs[1].plot(freqs * 1e-6, 20 * np.log10(np.abs(Pxx)))
    axs[1].set_xlabel('MHz')
    axs[1].set_title(r'$| P(\omega) |^2$')
    axs[1].set_ylabel('FFT Amplitude')

    fig.tight_layout()
    plt.show()


def compute_diff(t, pulse, opts, tr):
    """Compute the difference"""

    dx = t[1] - t[0]
    pdiff = np.gradient(pulse, dx)
    pd2 = pdiff ** 2

    G4 = np.trapz(pd2, dx=dx)

    print(f"G4 = {G4:.3e}")
    print(f"T_r = {tr:.3e}")
    print(f"T_r * G4 = {tr * G4:.2f}")

    #  plt.semilogy(t, pulse + 0.001, t, pdiff**2 + 0.001)
    #  plt.plot(t, pulse, t, pdiff**2)
    plt.plot(t, pulse / np.max(pulse), t, pd2 / np.max(pd2))
    plt.legend(['Pulse', r'$\dot{p}(t)^2_\mathrm{norm}$'])
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Look at some CRLB for TOA',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--sample_rate', action='store', type=float, default=200e6,
                        help='Sample rate')
    parser.add_argument('-d', '--duration', action='store', type=float, default=1e-6,
                        help='Duration of pulse')
    parser.add_argument('-t', '--type', action='store', choices=['lin', 'sin'], default='sin',
                        help='Type of rising edge (sin uses modified portions of rising/falling of sin)')
    parser.add_argument('-r', '--rise_percent', action='store', type=float, default=15,
                        help='Percent of the pulse duration to use as the rise time')
    args = parser.parse_args()

    pprint(args)
    t, p = main(args)

    #  plot_freq(t, p, args)
