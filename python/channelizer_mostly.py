

def _PrepareData(self, data):
    """docstring"""

    half_channels = self.num_channels >> 1
    new_len = int(int(data.size / half_channels) * half_channels) + 1
    new_len = new_len - half_channels if new_len > len(data) else new_len
    N = new_len
    num_columns = np.ceil((N - 1) / half_channels).astype(int) + 2

    D = np.zeros((self.num_channels, num_columns))
    D[:half_channels, 1:-1] = np.flipud(np.reshape(data[1:N], (half_channels, num_columns - 2), order='F'))
    D[half_channels:, 2:] = np.flipud(np.reshape(data[1:N], (half_channels, num_columns - 2), order='F'))
    D[0, 0] = data[0]
    D[half_channels, 1] = data[0]

    return D


def Channelize(self, data):
    """docstring"""

    # This takes the 1d input and does a couple of reshapes
    # and makes it into something that's *about* the size of:
    # (num_channels x data_len/decimation_factor). It's not incredibly cheap as is, maybe it
    # can be improved and help. But this function (Channelize) is about 17x more expensive
    # and so the place to start.
    poly_data = self._PrepareData(data)

    t = np.empty_like(poly_data)
    array_len = t.shape[1]
    num_channels = t.shape[0]
    for ch in range(num_channels):
        t[ch, :] = signal.convolve(poly_data[ch, :], self.poly_filt[ch, :], mode='full')[:array_len]

    y = np.fft.ifft(t, axis=0)

    # Perform modulation required on every other channel by fs/2
    v = y[:]
    v[1::2, 1::2] = -v[1::2, 1::2]

    v *= num_channels

    # Get rid of group delay (plus half a sample I've accounted for in xstart after the fact)
    v = v[:, self.out_skip:]

    return v


