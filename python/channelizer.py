#!/usr/bin/env python3
"""This is both a class and a test script for the Dave Sharpin version of the 2x
oversampled channelizer.

Trying to adjust this to be the Crochiere and Rabiner version instead of
Sharpin. Because the C&R one requires no shuffling of input data above what
can be done via strides and distances in FFTW.

For M = 3, K = 6, I = 2
Sharpin Forumalation:
    Data & Filter:
        x_p(n) = x(nM - p)
        h(nM + p) for n=2I, else 0

        x0  x3  x6  x9  x12 x15  |  h0  0   h6  0   h12 0
        x-1 x2  x5  x8  x11 x14  |  h1  0   h7  0   h13 0
        x-2 x1  x4  x7  x10 x13  |  h2  0   h8  0   h14 0
        x-3 x0  x3  x6  x9  x12  |  h3  0   h9  0   h15 0
        x-4 x-1 x2  x5  x8  x11  |  h4  0   h10 0   h16 0
        x-5 x-2 x1  x4  x7  x10  |  h5  0   h11 0   h17 0

Crochiere & Rabiner Forumulation:
    Data:
        x_p(n) = x(nM + p)
        h(nM - p) for n=2I, else 0

        x0  x3  x6  x9  x12 x15  |  h0  0   h6  0   h12 0
        x1  x4  x7  x10 x13 x16  |  h-1 0   h5  0   h11 0
        x2  x5  x8  x11 x14 x17  |  h-2 0   h4  0   h10 0
        x3  x6  x9  x12 x15 x18  |  h-3 0   h3  0   h9  0
        x4  x7  x10 x13 x16 x19  |  h-4 0   h2  0   h8  0
        x5  x8  x11 x14 x17 x20  |  h-5 0   h1  0   h7  0

"""

from copy import deepcopy
import numpy as np
from scipy import signal
from scipy.fft import rfft
import matplotlib.pyplot as plt

# {{{ class Channelizer
class Channelizer:
    """Class for the 2x oversampled channelizer structure in the Sharpin paper."""

    # {{{ __init__
    def __init__(self, sample_rate, decimation, overlap):
        self.sample_rate = sample_rate
        self.decimation = decimation
        self.overlap = overlap
        self.passband_ripple_db = 0.1
        self.stopband_rejection_db = -60
        self.out_bw = sample_rate / decimation / 2
        self.num_channels = decimation * 2
        self._design_filter()
        delta = self.sample_rate / self.decimation / 2

        self.freq_list = [np.arange(delta, self.sample_rate / 2, delta)]

        self.proto_filt_len = None
        self.proto_file = None
        self.grp_delay = None
        self.out_skip = None
        self.delay_delta = None

    # }}}

    # {{{ _design_filter
    def _design_filter(self):
        stop_edge = self.sample_rate / (2 * self.decimation)
        pass_edge = self.sample_rate / (4 * self.decimation) + self.overlap / 2
        tmp_fs = self.sample_rate

        if pass_edge >= stop_edge:
            raise RuntimeError(f"Passband edge {pass_edge} can't be higher than stop edge {stop_edge}")

        delta_f = stop_edge - pass_edge

        pass_ripple = min(
                [
                    10 ** (self.passband_ripple_db / 20) - 1,
                    1 - 10 ** (-self.passband_ripple_db / 20),
                ]
        )
        stop_ripple = 10 ** (self.stopband_rejection_db / 20)
        filt_len = (-2 / 3 * np.log10(10 * pass_ripple * stop_ripple) * tmp_fs / delta_f) / self.num_channels
        ideal_len = int(np.ceil(filt_len) * self.num_channels)
        self.proto_filt_len = int(np.floor(filt_len) * self.num_channels) + 1
        assert self.proto_filt_len <= ideal_len, "Prototype filter ideal length should be less than total length"

        amps = [1, 0]
        ripple = [1 / pass_ripple, 1 / stop_ripple]
        band_edges = [0, pass_edge, stop_edge, tmp_fs / 2]
        self.proto_filt = np.zeros(ideal_len)
        proto_filt = signal.remez(self.proto_filt_len, band_edges, amps, ripple, fs=tmp_fs)
        self.proto_filt[:len(proto_filt)] = proto_filt
        print(f"Prototype filter non-zero length: {self.proto_filt_len}, Padded length: {ideal_len}")

        #  self.grp_delay = (self.proto_filt.size - 1) / (2 * tmp_fs)
        self.grp_delay = (self.proto_filt_len - 1) / (2 * tmp_fs)
        self.out_skip = np.ceil(self.grp_delay * tmp_fs / self.decimation).astype(int)
        self.delay_delta = self.out_skip / (tmp_fs / self.decimation) - self.grp_delay
        print(f"Group delay: {self.grp_delay}, output samples: {self.grp_delay * tmp_fs / self.decimation}")

        up_size = int(len(self.proto_filt) / self.num_channels * 2)
        self.poly_filt = np.zeros((self.num_channels, up_size))

        decimation = self.num_channels >> 1
        for idx in range(self.poly_filt.shape[1]):
            for chan in range(self.poly_filt.shape[0]):
                in_idx = idx * decimation - chan
                if idx % 2 == 0 and 0 <= in_idx < len(self.proto_filt):
                        self.poly_filt[chan, idx] = self.proto_filt[in_idx]


    # }}}

    # {{{ _prepare_data
    def _prepare_data(self, data):
        """Ignore last few samples"""

        decimation = self.num_channels >> 1
        new_len = (len(data) // decimation) * decimation
        num_columns = ((new_len - 1) // decimation) + 1
        data_poly = np.zeros((self.num_channels, num_columns))
        tmp_data = np.reshape(data[:new_len], (decimation, -1), order="F")
        data_poly[:decimation, :] = tmp_data
        data_poly[decimation:, :-1] = deepcopy(tmp_data[:, 1:])

        return data_poly

    # }}}

    # {{{ channelize
    def channelize(self, data):
        """Perform the actual channelization step"""

        poly_data = self._prepare_data(data)

        filt_out = np.empty_like(poly_data)
        array_len = filt_out.shape[1]
        for chan in range(self.num_channels):
            filt_out[chan, :] = signal.convolve(poly_data[chan, :], self.poly_filt[chan, :], mode="full")[:array_len]

        fft_out = np.conj(rfft(filt_out, axis=0, overwrite_x=True))

        mod_out = fft_out[:]
        mod_out[1::2, 1::2] = -mod_out[1::2, 1::2]

        data_channelized = mod_out[:, self.out_skip:]

        return data_channelized[1 : self.decimation]

    # }}}

# }}}

# {{{ __main__
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Test Channelizer", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sample-rate", help="Sample rate (Hz)", action="store", type=float, dest="sample_rate", default=800)
    parser.add_argument("--overlap", help="Overlap (Hz)", action="store", type=float, default=1)
    parser.add_argument("--decimation", help="Decimation factor", action="store", type=int, default=32)
    parser.add_argument("--num-points", help="Number of test input points", action="store", type=int, dest="num_points", default=100000)
    opts = parser.parse_args()

    time = np.arange(opts.num_points) / opts.sample_rate
    num_cycles = 2
    sig = signal.chirp(time, 0, time[-1], opts.sample_rate * num_cycles)

    chann = Channelizer(opts.sample_rate, opts.decimation, opts.overlap)

    v_sig = chann.channelize(sig)

    chan_max = v_sig.shape[0]
    time_max = v_sig.shape[1]

    out_per = opts.decimation / opts.sample_rate
    time_series = np.arange(time_max + 1) * out_per + out_per / 2
    channs = np.arange(chan_max + 1) + 0.5

    plt.rcParams["axes.grid"] = False
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].grid(False)
    axs[0, 0].pcolormesh(time_series, channs, np.abs(v_sig))
    axs[0, 0].set_xlabel("Time")
    axs[0, 0].set_ylabel("Channel")
    axs[0, 0].set_title("Channelizer Output Power")

    axs[1, 0].grid(True)
    d = v_sig[1:3, :]
    dd = np.angle(d[:, 1:] * np.conj(d[:, :-1]))
    axs[1, 0].plot(time_series[:-2], dd.T)
    axs[1, 0].set_xlabel("Time")
    axs[1, 0].set_ylabel("Channel")
    axs[1, 0].set_title("Two Channels Phase DCM")
    axs[1, 0].grid(True)

    axs[0, 1].plot(time_series[:-1], 20 * np.log10(np.abs(v_sig).T))
    axs[0, 1].set_title("Amplitude per Channel")
    axs[0, 1].set_xlabel("Time")
    axs[0, 1].set_ylabel("Amplitude (dB)")
    axs[0, 1].grid(True)

    axs[1, 1].plot(chann.poly_filt.T)
    axs[1, 1].set_title("Polyphase Filter Decomp")
    axs[1, 1].set_xlabel("Time")
    axs[1, 1].set_ylabel("Amplitude")
    axs[1, 1].grid(True)
    axs[1, 1].legend()

    plt.tight_layout()
    plt.show()

# }}}
