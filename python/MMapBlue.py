#!/usr/bin/env python3

import numpy as np

class MMapBlue:
    def __init__(self, filename, ext_header_size):
        '''Memory map input file and set up limits on reading'''
        self.ext_header_size = ext_header_size
        self.a = np.arange(131072 + 512 + self.ext_header_size, dtype=np.int8)

    def __len__(self):
        '''Length of the array'''
        return len(self.a) - 512 - self.ext_header_size

    def __getitem__(self, idx):
        '''Get index'''
        if idx >= len(self.a) - 512 - self.ext_header_size:
            raise IndexError(f"Index of {idx} is out of bounds for axis 0")

        return self.a[idx - 512]


