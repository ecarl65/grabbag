#!/usr/bin/env python3

import numpy as np
import DirectFiltering
import matplotlib.pyplot as plt
import scipy as sp
import scipy.signal as signal


if __name__ == "__main__":
    h0 = signal.firwin(33, 0.2, pass_zero=True)
    h0z = np.zeros(2 * len(h0) - 1)
    h0z[::2] = h0
    s = np.random.randn(1000) * 0.1 + np.sin(np.arange(1000) * 0.1)

    c = DirectFiltering.DirectFIRZeros(s, h0z)

    t = np.arange(len(s))

    plt.plot(t, s, label='Original')
    plt.plot(t, c, label='Upzero filtered')
    plt.legend()
    plt.show()
