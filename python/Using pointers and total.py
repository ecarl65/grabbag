    cdef double *pdata = <double*>&data[0]
    cdef double *pfir = <double*>&fir[0]
    cdef double *presult = <double*>&result_view[0]
    cdef double total = 0.0

    # Do filtering in three sections. Boundary conditions on beginning and end
    # with the bulk being the middle section without boundary conditions
    for n in range (half_len):
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
