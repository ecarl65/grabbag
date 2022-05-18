#!/usr/bin/env python3

from copy import deepcopy
import numpy as np
from scipy import signal
from scipy.fft import rfft
import matplotlib.pyplot as plt

# {{{ class Channelizer
class Channelizer:
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
        ideal_size = 2 * self.decimation
        filt_len = (-2 / 3 * np.log10(10 * pass_ripple * stop_ripple) * tmp_fs / delta_f) / ideal_size
        self.proto_filt_len = int(np.ceil(filt_len) * ideal_size)

        amps = [1, 0]
        ripple = [1 / pass_ripple, 1 / stop_ripple]
        band_edges = [0, pass_edge, stop_edge, tmp_fs / 2]
        self.proto_filt = signal.remez(self.proto_filt_len, band_edges, amps, ripple, fs=tmp_fs)

        self.grp_delay = (self.proto_filt.size - 1) / (2 * tmp_fs)
        self.out_skip = np.ceil(self.grp_delay * tmp_fs / self.decimation).astype(int)
        self.delay_delta = self.out_skip / (tmp_fs / self.decimation) - self.grp_delay

        e_cs = np.reshape(
                self.proto_filt,
                (self.num_channels, int(len(self.proto_filt) / self.num_channels)),
                order="F",
        )
        up_size = int(len(self.proto_filt) / self.num_channels * 2)
        self.poly_filt = np.zeros((self.num_channels, up_size))
        self.poly_filt[:, ::2] = e_cs

    # }}}

    # {{{ _prepare_data
    def _prepare_data(self, data):

        half_channels = self.num_channels >> 1
        new_len = int(int(data.size / half_channels) * half_channels) + 1
        new_len = new_len - half_channels if new_len > len(data) else new_len
        num_columns = np.ceil((new_len - 1) / half_channels).astype(int) + 2

        data_poly = np.zeros((self.num_channels, num_columns))
        tmp_data = np.flipud(np.reshape(data[1:new_len], (half_channels, num_columns - 2), order="F"))
        data_poly[:half_channels, 1:-1] = tmp_data
        data_poly[half_channels:, 2:] = deepcopy(tmp_data)
        data_poly[0, 0] = data[0]
        data_poly[half_channels, 1] = data[0]

        return data_poly

    # }}}

    # {{{ channelize
    def channelize(self, data):

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

    # {{{ channelize_chunks
    def channelize_chunks(self, data, outfile):
        """Break channelizer into chunks of data and write output to file"""

        if len(data) < 4 * self.proto_filt_len:
            outdata = self.channelize(data)
            # TODO Write outdata to outfile all at once
            return True

        # Since returned above don't need "else", reduce branch depth
        num_overlap_save = 4 * self.proto_filt_len
        num_valid = num_overlap_save - self.proto_filt_len + 1 # XXX Verify +1
        num_buffers = len(data) // num_valid

        for buff in range(num_buffers):
            strt = num_valid * buff
            stp = strt + num_overlap_save

            data_chunk = data[strt:stp]

            chunk_out = self.channelize(data_chunk)

            # Figure out what to save
            # If it were just FFT stuff it'd be
            # chunk_out[self.proto_filt_len - 1 : ]
            # But I don't think that's right for the channelizer



        # Iterate through chunks and do some overlap

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

    # plot output
    #  plt.plot(np.abs(v_sig.T))

    chan_max = v_sig.shape[0]
    time_max = v_sig.shape[1]

    delta = opts.decimation / opts.sample_rate
    time_series = np.arange(time_max + 1) * delta + delta / 2
    channs = np.arange(chan_max + 1) + 0.5

    plt.pcolormesh(time_series, channs, np.abs(v_sig))
    plt.xlabel("Time")
    plt.ylabel("Channel")
    plt.title("Channelizer Output Power")

    plt.show()

# }}}