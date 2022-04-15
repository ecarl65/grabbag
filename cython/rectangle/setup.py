from distutils.core import setup, Extension
import numpy as np
from Cython.Build import cythonize

setup(
        ext_modules=cythonize("rect.pyx"),
        include_dirs=[np.get_include()]
        )

#  setup(ext_modules=[
    #  Extension("rect.pyx", ["rect.c"],
        #  include_dirs=[np.get_include()]),
    #  ],
#  )
