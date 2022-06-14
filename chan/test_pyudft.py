#!/usr/bin/env python3
### current bpython session - make changes and save to reevaluate session.
### lines beginning with ### will be ignored.
### To return to bpython without reevaluating make no changes to this file
### or save an empty file.
import numpy as np
import pyudft
from scipy import signal
import matplotlib.pyplot as plt

M = 8
Nfilt = M * 8 + 1
Fs = 10e3
write_out = True
verbose = True
pu = pyudft.pudft(M, Nfilt, Fs, write_out, verbose)
t = np.arange(2048) / Fs
b = signal.chirp(t, 0, t[-1], 2 * Fs).astype('f')
od = pu.pyrun(b)
od
# OUT: array([[ 0.00000000e+00+0.00000000e+00j,  0.00000000e+00+0.00000000e+00j,
# OUT:          0.00000000e+00+0.00000000e+00j,  0.00000000e+00+0.00000000e+00j,
# OUT:          0.00000000e+00+0.00000000e+00j],
# OUT:        [ 0.00000000e+00+0.00000000e+00j,  0.00000000e+00+0.00000000e+00j,
# OUT:          0.00000000e+00+0.00000000e+00j,  0.00000000e+00+0.00000000e+00j,
# OUT:          0.00000000e+00+0.00000000e+00j],
# OUT:        [ 0.00000000e+00+0.00000000e+00j,  0.00000000e+00+0.00000000e+00j,
# OUT:          0.00000000e+00+0.00000000e+00j,  0.00000000e+00+0.00000000e+00j,
# OUT:          0.00000000e+00+0.00000000e+00j],
# OUT:        ...,
# OUT:        [-4.08494882e-02+0.00000000e+00j, -4.71817236e-03+6.58320729e-03j,
# OUT:         -7.69061968e-04-1.06798485e-04j, -4.40319302e-04-2.49922741e-05j,
# OUT:         -3.73832881e-04+0.00000000e+00j],
# OUT:        [ 7.52543807e-01+0.00000000e+00j,  6.06995821e-03+1.12580357e-03j,
# OUT:          2.31546164e-03+2.99662352e-04j,  1.32966018e-03+6.85959822e-05j,
# OUT:          1.12995505e-03+0.00000000e+00j],
# OUT:        [ 1.06223488e+00+0.00000000e+00j, -2.28903145e-02-1.21171894e-02j,
# OUT:         -6.79385662e-03-1.11514330e-03j, -3.94779351e-03-2.62183137e-04j,
# OUT:         -3.36390734e-03+0.00000000e+00j]])
plt.plot(np.abs(od)); plt.show()
# OUT: [<matplotlib.lines.Line2D object at 0x7f37a3fba320>, <matplotlib.lines.Line2D object at 0x7f37a3fba2c0>, <matplotlib.lines.Line2D object at 0x7f37a3fba3e0>, <matplotlib.lines.Line2D object at 0x7f37a3fba500>, <matplotlib.lines.Line2D object at 0x7f37a3fba620>]
### 
