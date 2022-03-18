import math
import numpy

def zero_padded_grab(data, offset, length):
    """This helper function allows grabbing (copying or creating a view)
    of data before, during and after the array.  Out of bounds samples
    are filled with zeros.
    """

    assert length > 0, "length must be non-negative"

    # fast case first, return low overhead view
    if offset >= 0 and offset + length <= len(data):
        return data[offset : offset + length]

    # slow copy and concat
    if offset < 0:
        amount = min(-offset, length)
        before = numpy.zeros(amount)
        length -= amount
        offset += amount
    else:
        before = numpy.zeros(0)

    if offset < len(data):
        amount = min(length, len(data) - offset)
        during = data[offset : offset + amount]
        length -= amount
    else:
        during = numpy.zeros(0)

    return numpy.concatenate((before, during, numpy.zeros(length)))

def realchan(data, fir, chans):
    """Performs a polyphase channelizer on the real-valued input data
    using the provided fir as a filter.  The data and fir arrays must
    have one dimension, and the return value will be a complex-valued
    two dimensional array of data where the columns are the separate
    channels, and the rows increase in the "time" direction.

    Each column is a separate frequency using a real FFT
    (see numpy.fft.rfft), and this will include DC and Nyquist.

    As an example, if you request 5 channels, you will get the
    following frequencies:  [DC, 1*Fs/8, 2*Fs/8, 3*Fs/8, Nyquist]
    """

    # These won't make copies if the arguments are already arrays
    data = numpy.asarray(data)
    fir = numpy.asarray(fir)

    assert len(data.shape) == 1, "expecting 1D data"
    assert len(fir.shape) == 1, "expecting 1D fir filter"
    assert len(fir)%2 == 1, "expecting an odd-length fir filter"
    assert chans >= 2, "expecting at least 2 channels"
    assert numpy.isrealobj(data), "expecting real data"
    assert numpy.isrealobj(fir), "expecting real fir filter"

    cols = (chans - 1)*2
    half_length = len(fir)//2
    rows_before = math.ceil(half_length/cols)
    rows_after = math.ceil((half_length + 1 + cols//2)/cols)
    filter_rows = rows_before + rows_after

    # We're doing FFT based convolution vertically, so we want
    # to pick an FFT size big enough to be worth the overhead
    rows = 32
    while rows < 4*filter_rows: rows *= 2

    filter_offset = half_length - rows_before*cols
    even_filter = zero_padded_grab(fir, filter_offset, rows*cols)
    odd_filter = zero_padded_grab(fir, filter_offset - cols//2, rows*cols)
    even_filter = even_filter.reshape((rows, cols))
    odd_filter = odd_filter.reshape((rows, cols))

    # do vertical FFTs down the filters and scale
    filter_scale = 1.0/(rows//2 + 1)
    even_filter = numpy.fft.rfft(even_filter, axis=0).conj()*filter_scale
    odd_filter = numpy.fft.rfft(odd_filter, axis=0).conj()*filter_scale

    data_offset = -rows_before*cols
    row_progress = rows - filter_rows + 1
    data_rows = math.ceil(len(data)/cols)
    iterations = math.ceil(data_rows/row_progress)
    result = numpy.empty(
        shape=(iterations*2*row_progress, chans),
        dtype=numpy.complex128
    )

    for ii in range(iterations):
        grab = zero_padded_grab(data, data_offset, rows*cols)
        grab = grab.reshape((rows, cols))
        grab = numpy.fft.rfft(grab, axis=0)
        even = numpy.fft.irfft(grab*even_filter, axis=0)
        odd = numpy.fft.irfft(grab*odd_filter, axis=0)
        even = numpy.fft.rfft(even, axis=1)
        odd = numpy.fft.rfft(odd, axis=1)
        start = ii*row_progress*2
        stop = start + 2*row_progress
        result[start + 0 : stop + 0 : 2, :] = even[:row_progress,:]
        result[start + 1 : stop + 1 : 2, :] = odd[:row_progress,:]
        data_offset += row_progress*cols

    return result

#
# Visualization for chans == 4, len(fir) == 15
#
# data: 0 0 0 0 0 0 0 0 0 0 a b c d e f g h i j k l m n o p q r s t u v w x y z
#                           |
#                           |  <- this is the middle of the FIR
#                           |
# fir:        A B C D E F G H I J K L M N O
#
# Reshaped:
#         data:               even fir:          odd filter:
#            0 0 0 0 0 0         0 0 0 0 0 A        0 0 0 0 0 0
#            0 0 0 0 0 0         B C D E F G        0 0 A B C D
#            a b c d e f         H I J K L M        E F G H I J
#            g h i j k l         N O 0 0 0 0        K L M N O 0
#            m n o p q r         0 0 0 0 0 0        0 0 0 0 0 0
#            s t u v w x         0 0 0 0 0 0        0 0 0 0 0 0
#            y x . . . .         0 0 0 0 0 0        0 0 0 0 0 0
#            . . . . . .         0 0 0 0 0 0        0 0 0 0 0 0
#            . . . . . .         0 0 0 0 0 0        0 0 0 0 0 0
#            . . . . . .         0 0 0 0 0 0        0 0 0 0 0 0
#                                 


################ Everything below here is testing/demonstration

from math import cos, sin, pi
from os import popen
from time import time

def demo():
    bins = 64
    samplerate = 25e6
    bandwidth = samplerate/bins
    transition = .2*bandwidth
    seconds = 10.0
    channels = bins//2 + 1
    testfreq = .125*samplerate + 1
    testphase = 60.0 # degrees

    # build the data and fir
    length = math.ceil(samplerate*seconds)
    print("Building data...")
    ramp = numpy.arange(length)
    data = numpy.cos(2*pi*(testfreq/samplerate*ramp + testphase/360))
    #data = 0.0*ramp
    #data[len(ramp)//4] = 1.0
    print("Building fir...")
    fir = lowpass(2, samplerate, bandwidth, transition)
    print("len(fir) ==", len(fir))

    # run he channelizer and save as type 2000
    print("Running realchan...")
    before = time()
    output = realchan(data, fir, channels)
    print("channelizer took", time() - before, "seconds")
    print("output.shape:", output.shape)
    print("expected:", len(data)/(2 * (channels-1))*2)
    dump2000("channels.tmp", output)

    if 0:
        gnuplot = popen("gnuplot -persist", 'w')
        gnuplot.write("plot '-' w lines t 'fir'\n")
        for ii in range(len(fir)):
            gnuplot.write("%f\n" % fir[ii]);
        gnuplot.write("e\n")
        gnuplot.flush()

    if 1:
        padded = numpy.concatenate((fir, numpy.zeros(len(fir)*20)))
        spec = numpy.fft.rfft(padded)
        spec = numpy.abs(spec)
        spec = 20*numpy.log10(spec + 1e-10)
        gnuplot = popen("gnuplot -persist", 'w')
        gnuplot.write("plot '-' w lines t 'freq response'\n")
        length = len(spec)
        for ii in range(len(spec)):
            freq = .5*ii/length
            gnuplot.write("%f %f\n" % (freq, spec[ii] - spec[0]));
        gnuplot.write("e\n")
        gnuplot.flush()

def sinc(xx):
    if xx == 0.0: return 1.0
    return sin(pi*xx)/(pi*xx)

def firwin(window, offset, length):
    if offset > +.5*length: return 0;
    if offset < -.5*length: return 0;

    xx = 2.0*pi*offset/length
    if window == 0: return 1.0
    if window == 1: return 0.5 + 0.5*cos(xx)
    if window == 2: return .41 + 0.5*cos(xx) + 0.09*cos(2.0*xx)

    raise ValueError("unsupported window")

def firtaps(window, inrate, transbw):
    if window == 0: return 1.000000*inrate/transbw
    if window == 1: return 1.673904*inrate/transbw
    if window == 2: return 2.810852*inrate/transbw

    raise ValueError("unsupported window")

def lowpass(window, samplerate, bandwidth, transbw):
    taps = math.ceil(firtaps(1, samplerate, transbw))
    if taps%2 == 0: taps += 1
    half = taps//2
    width = bandwidth/samplerate
    fir = numpy.empty((taps,), dtype=numpy.float64)
    for ii in range(taps):
        xx = ii - half
        fir[ii] = sinc(xx*width)*firwin(window, xx, taps)
    return fir

def dump2000(path, data):
    assert data.dtype == numpy.complex128
    ycount, xcount = data.shape
    import struct
    xdelta = .5/(xcount - 1)
    print("xcount:", xcount, "xdelta:", xdelta)
    header = struct.pack(
        '<4s4s4s20xddi2s2xd192xddiiddii208x',
        b"BLUE", b"EEEI", b"EEEI", 512.0,
        xcount*ycount*16.0, 2000, b"CD", 0.0,
        0.0, xdelta, 0, xcount,
        0.0, 1.0, 0, 0
    )
    print("len(header):", len(header))
    with open(path, 'w+b') as f:
        f.write(header)
        f.write(data)

if __name__ == '__main__':
    demo()

    
