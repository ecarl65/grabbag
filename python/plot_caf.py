#!/usr/bin/env python3 

import matplotlib.pyplot as plt
import numpy as np
import bluefile as bf
import sys

infile = sys.argv[1]

bf.set_type2000_format(np.ndarray)
h, d = bf.read(infile)
fdoa = np.arange(len(d[0]) + 1) * h['xdelta'] + h['xstart']
tdoa = np.arange(len(d) + 1) * h['ydelta'] + h['ystart']

cmap = plt.colormaps['gray']
plt.grid(False)

#  plt.pcolormesh(fdoa, 1e3 * tdoa, np.abs(d), shading='flat', cmap=cmap)
plt.pcolor(fdoa, 1e3 * tdoa, np.abs(d), shading='flat', cmap=cmap)
plt.xlabel('FDOA (Hz)')
plt.ylabel('TDOA (ms)')

plt.show()
