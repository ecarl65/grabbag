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
    return 20 * np.log10(np.abs(x) + 1e-8)

def p2db(x):
    return 10 * np.log10(np.abs(x) + 1e-8)

# {{{ Noise
class Noise:

    # {{{ __init__
    def __init__(self, num_samps, samp_rate, snr, nfft, pad, no_plots=False, psd_type="psd"):
        self.num_samps = num_samps
        self.num_signal = num_samps // 2
        self.samp_rate = samp_rate
        self.snr = snr
        self.nfft = nfft
        self.overlap = nfft // 2
        self.pad = pad
        self.no_plots = no_plots
        self.psd_type = psd_type

        # Input signal and noise
        self.time = np.arange(self.num_samps) / self.samp_rate
        self.freq_0 = self.samp_rate / self.nfft * 10.00  # Alter fractional part to hit bin center or not
        self.complex_noise = cnormal(len(self.time)) * 10 ** (-self.snr / 20)
        self.complex_signal = np.exp(1j * 2 * np.pi * self.freq_0 * self.time)
        # self.complex_signal = 1
        if self.psd_type == "stft":
            self.complex_signal *= np.r_[np.ones(self.num_signal), np.zeros(self.num_samps - self.num_signal)]
        self.spn = self.complex_signal + self.complex_noise

        # Expected gains and SNRs
        self.proc_gain = p2db(self.nfft)
        self.expect = self.proc_gain + self.snr

        # Better to do both sides but for now just get noise stats on one side of signal
        self.ncutoff = int(self.pad * 0.4)

        self.nperseg = 16

    # }}}

    # {{{ run_psd
    def run_psd(self):

        # To normalize all outputs the same treat FFT like multiplying by rect window and normalize "filter"
        h_rect = np.ones(self.nfft) / self.nfft
        SPN = np.fft.fftshift(np.fft.fft(h_rect * self.spn[:self.nfft], self.pad))

        # My own periodogram. Normalize filter as well for same peak power. No overlap.
        win = signal.windows.hamming(self.nfft)
        win /= np.sum(win)
        n_loops = int(self.num_samps / self.overlap) - 1
        SPN_per = np.zeros(self.pad)
        new_samps = self.nfft - self.overlap
        for iter in range(n_loops):
            data = self.spn[iter * new_samps : iter * new_samps + self.nfft]
            SPN_per += np.abs(np.fft.fftshift(np.fft.fft(win * data, self.pad)))**2
        SPN_per /= n_loops

        # Matlab PSD
        Pxx, freqs = psd(self.spn, Fs=self.samp_rate, NFFT=self.nfft, sides="twosided", noverlap=self.overlap,
                         pad_to=self.pad, scale_by_freq=False)

        # Plot results
        if not self.no_plots:
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

        # Print out some results
        print(f"Input SNR: {self.snr}")
        print(f"Processing gain: {self.proc_gain:.2f}")
        print(f"Expected SNR: {self.expect:.2f}")
        print(f"Number of averages in periodogram: {n_loops}")

        # Quantify SNR per method
        noise_mean = {
            "PSD": np.mean(Pxx[:self.ncutoff]),
            "FFT": np.mean(np.abs(SPN[:self.ncutoff])**2),
            "Per": np.mean(SPN_per[:self.ncutoff]),
            #  "PSD": np.median(Pxx) / np.log(2),
            #  "FFT": np.median(np.abs(SPN)**2) / np.log(2),
            #  "Per": np.median(SPN_per) / np.log(2),
        }
        for psd_type, nm in noise_mean.items():
            snr_db = p2db((1 - nm) / nm)
            # snr_db = p2db(1 / nm)
            print(f"{psd_type}: SNR = {snr_db:.2f}")
            print(f"{psd_type}: Delta = {self.expect - snr_db:.2f}")

    # }}}

    # {{{ run
    def run(self):
        if self.psd_type == "psd":
            self.run_psd()
        elif self.psd_type == "stft":
            self.run_stft()
        else:
            raise RuntimeError("No valid run method chosen")
    # }}}

    # {{{ run_stft
    def run_stft(self):

        # Compute STFT
        #  win = np.ones(self.nperseg) / self.nperseg
        win = signal.windows.hann(self.nperseg)
        win /= np.sum(win)
        freqs, times, stft = signal.stft(self.spn, window=win, fs=self.samp_rate, nperseg=self.nperseg, return_onesided=False)
        #  freqs, times, stft = signal.stft(self.spn, window="hamming", fs=self.samp_rate, nperseg=self.nperseg, return_onesided=False)
        stft[:, 1::2] *= -1
        stft = np.fft.fftshift(stft, axes=0)
        freqs = np.fft.fftshift(freqs)

        # Compute stats
        full_var = np.var(self.complex_noise)
        out_var  = np.var(stft[:, stft.shape[1]//2 + 1:], axis=1)
        filt_freqs, freq_response = signal.freqz(win, 1, fs=self.samp_rate)
        neb = np.trapz(np.abs(freq_response / np.max(freq_response)) ** 2, filt_freqs)
        print(f"Ratio of input noise std to outout: {full_var / np.median(out_var)}")
        print(f"Ratio of half sample rate to NEB cutoff: {self.samp_rate / 2 / neb}")

        # Plot the data in various ways
        if not self.no_plots:
            plt.close("all")
            fig, axs = plt.subplots(2, figsize=(24,12), tight_layout=True)
            axs[0].plot(self.time, a2db(self.complex_signal), color="k", linewidth=2, label="Full Rate Signal")
            axs[0].plot(self.time, a2db(self.complex_noise), color="dimgray", linewidth=2, label="Full Rate Noise", alpha=0.5)
            axs[0].plot(times, a2db(stft.T), alpha=1.0)
            axs[0].set_title("Time Domain STFT Output")
            axs[0].set_xlabel("Time (sec)")
            axs[0].set_ylabel("Power (dB)")
            axs[0].legend()

            axs[1].grid(False)
            pcm = axs[1].pcolormesh(times, freqs, a2db(stft), shading="auto")
            axs[1].set_xlabel("Time (sec)")
            axs[1].set_ylabel("Freq (Hz)")
            axs[1].set_title("STFT Power (dB)")
            fig.colorbar(pcm, ax=axs[1])
            plt.show()

    # }}}

# }}}

# {{{ __main__
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Test SNR", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sample-rate", help="Sample rate (Hz)", action="store", type=float, dest="sample_rate", default=2.5e6)
    parser.add_argument("--num-samples", help="Number of samples", action="store", type=float, default=8192)
    parser.add_argument("--nfft", help="Number of input data samples per FFT", type=int, default=1024)
    parser.add_argument("-p", "--pad", help="Size of FFT, greater than or equal to nfft", type=int, default=1024)
    parser.add_argument("-s", "--snr", help="SNR in dB", type=float, default=20)
    parser.add_argument("-n", "--no-plot", help="Set to suppress plotting", action="store_true")
    parser.add_argument("-t", "--type", help="Which method to run", choices=["psd", "stft"], default="psd")
    opts = parser.parse_args()

    sig_noise = Noise(opts.num_samples, opts.sample_rate, opts.snr, opts.nfft, opts.pad, opts.no_plot, opts.type)
    sig_noise.run()

# }}}
