# cython: language_level=3, boundscheck=False, wraparound=False

import numpy as np

# {{{ DirectSymmetricFIR
def DirectSymmetricFIR(double[:] data, double[:] fir):
    """Perform a filtering operation on data with a symmetric, even length FIR filter

    Arguments:
        data - Input numpy array of doubles
        fir - Input numpy array of doubles that's an even length symmetric FIR filter

    Returns:
        result - Result of filtering"""

    cdef Py_ssize_t data_len = data.size
    cdef Py_ssize_t filt_len = fir.size
    cdef Py_ssize_t half_len = filt_len >> 1
    result = np.zeros(data_len)
    cdef double[:] result_view = result
    cdef double tmp_val = 0
    cdef Py_ssize_t tmp_idx = 0
    cdef Py_ssize_t n
    cdef Py_ssize_t i

    cdef double *pdata = <double*>&data[0]
    cdef double *pfir = <double*>&fir[0]
    cdef double *presult = <double*>&result_view[0]
    cdef double total = 0.0

    # Do filtering in three sections. Boundary conditions on beginning and end
    # with the bulk being the middle section without boundary conditions
    for n in range(half_len):
        total = 0.0
        for i in range(0, half_len):
            tmp_idx = n - half_len + i
            if tmp_idx < 0:
                tmp_val = 0
            else:
                tmp_val = pdata[tmp_idx]
            total += pfir[i] * (tmp_val + pdata[n + half_len - i + 1])
        presult[n] = total

    for n in range(half_len, data_len - half_len + 1):
        total = 0.0
        for i in range(0, half_len):
            total += pfir[i] * (pdata[n - half_len + i] + pdata[n + half_len - i + 1])
        presult[n] = total

    for n in range(data_len - half_len + 1, data_len):
        total = 0.0
        for i in range(0, half_len):
            tmp_idx = n + half_len - i + 1
            if tmp_idx >= data_len:
                tmp_val = 0
            else:
                tmp_val = pdata[tmp_idx]
            total += pfir[i] * (pdata[n - half_len + i] + tmp_val)
        presult[n]= total

    return result
# }}}


# {{{ DirectFIR
def DirectFIR(double[:] data, double[:] fir):
    """Perform a filtering operation on data with FIR filter of *odd* length

    Arguments:
        data - Input numpy array of doubles
        fir - Input numpy array of doubles that's an even length symmetric FIR filter

    Returns:
        result - Result of filtering"""

    cdef Py_ssize_t data_len = data.size
    cdef Py_ssize_t filt_len = fir.size
    cdef Py_ssize_t half_len = filt_len // 2
    if filt_len % 2 == 0:
        raise NotImplemented("This only handles odd length filters for now")
    result = np.zeros(data_len)
    cdef double[:] result_view = result
    cdef double tmp_val = 0
    cdef Py_ssize_t tmp_idx = 0
    cdef Py_ssize_t n
    cdef Py_ssize_t i

    cdef double *pdata = <double*>&data[0]
    cdef double *pfir = <double*>&fir[0]
    cdef double *presult = <double*>&result_view[0]
    cdef double total = 0.0

    # Do filtering in three sections. Boundary conditions on beginning and end
    # with the bulk being the middle section without boundary conditions
    for n in range(half_len):
        total = 0.0
        for i in range(filt_len):
            tmp_idx = n - i + half_len
            if tmp_idx < 0:
                continue
            total += pfir[i] * pdata[tmp_idx]
        presult[n] = total

    for n in range(half_len, data_len - half_len):
        total = 0.0
        for i in range(filt_len):
            total += pfir[i] * pdata[n - i + half_len]
        presult[n] = total

    for n in range(data_len - half_len, data_len):
        total = 0.0
        for i in range(filt_len):
            tmp_idx = n - i + half_len
            if tmp_idx >= data_len:
                tmp_val = 0
            else:
                tmp_val = pdata[tmp_idx]
            total += pfir[i] * tmp_val
        presult[n]= total

    return result
# }}}


# {{{ DirectFIRCopy
def DirectFIRCopy(double[:] data, double[:] fir, double[:] result):
    """Perform a filtering operation on data with FIR filter of *odd* length

    Arguments:
        data - Input numpy array of doubles
        fir - Input numpy array of doubles that's an even length symmetric FIR filter
        result - Output array to write to"""


    cdef Py_ssize_t out_len = result.size
    cdef Py_ssize_t data_len = data.size
    cdef Py_ssize_t filt_len = fir.size
    cdef Py_ssize_t half_len = filt_len // 2
    if filt_len % 2 == 0:
        raise NotImplemented("This only handles odd length filters for now")
    if data_len != out_len:
        raise Exception("Input array and output array should be same size")
    cdef double tmp_val = 0
    cdef Py_ssize_t tmp_idx = 0
    cdef Py_ssize_t n
    cdef Py_ssize_t i

    cdef double *pdata = <double*>&data[0]
    cdef double *pfir = <double*>&fir[0]
    cdef double *presult = <double*>&result[0]
    cdef double total = 0.0

    # Do filtering in three sections. Boundary conditions on beginning and end
    # with the bulk being the middle section without boundary conditions
    for n in range(half_len):
        total = 0.0
        for i in range(filt_len):
            tmp_idx = n - i + half_len
            if tmp_idx < 0:
                continue
            total += pfir[i] * pdata[tmp_idx]
        presult[n] = total

    for n in range(half_len, data_len - half_len):
        total = 0.0
        for i in range(filt_len):
            total += pfir[i] * pdata[n - i + half_len]
        presult[n] = total

    for n in range(data_len - half_len, data_len):
        total = 0.0
        for i in range(filt_len):
            tmp_idx = n - i + half_len
            if tmp_idx >= data_len:
                tmp_val = 0
            else:
                tmp_val = pdata[tmp_idx]
            total += pfir[i] * tmp_val
        presult[n]= total

    return result
# }}}


# {{{ DirectFIRZeros
def DirectFIRZeros(double[:] data, double[:] fir):
    """Perform a filtering operation on data with FIR filter of *odd* length and
    every other sample being zero

    Arguments:
        data - Input numpy array of doubles
        fir - Input numpy array of doubles that's an even length symmetric FIR filter

    Returns:
        result - Result of filtering"""

    cdef offset = 0
    cdef Py_ssize_t data_len = data.size
    cdef Py_ssize_t filt_len = fir.size
    cdef Py_ssize_t half_len = filt_len // 2
    if filt_len % 2 == 0:
        raise NotImplemented("This only handles odd length filters for now")
    result = np.zeros(data_len)
    cdef double[:] result_view = result
    cdef double tmp_val = 0
    cdef Py_ssize_t tmp_idx = 0
    cdef Py_ssize_t n
    cdef Py_ssize_t i

    cdef double *pdata = <double*>&data[0]
    cdef double *pfir = <double*>&fir[0]
    cdef double *presult = <double*>&result_view[0]
    cdef double total = 0.0

    # Check which sample the zeros start on
    if pfir[0] == 0.0:
        offset = 1
    elif pfir[1] == 0.0:
        offset = 0
    else:
        raise Exception("This function only works with filter that has every other sample as zero, starting with either first or second sample")

    # Do filtering in three sections. Boundary conditions on beginning and end
    # with the bulk being the middle section without boundary conditions
    for n in range(half_len):
        total = 0.0
        for i in range(offset, filt_len, 2):
            tmp_idx = n - i + half_len
            if tmp_idx < 0:
                continue
            total += pfir[i] * pdata[tmp_idx]
        presult[n] = total

    for n in range(half_len, data_len - half_len):
        total = 0.0
        for i in range(offset, filt_len, 2):
            total += pfir[i] * pdata[n - i + half_len]
        presult[n] = total

    for n in range(data_len - half_len, data_len):
        total = 0.0
        for i in range(offset, filt_len, 2):
            tmp_idx = n - i + half_len
            if tmp_idx >= data_len:
                tmp_val = 0
            else:
                tmp_val = pdata[tmp_idx]
            total += pfir[i] * tmp_val
        presult[n]= total

    return result
# }}}

