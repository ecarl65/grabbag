#!/usr/bin/env python

from setuptools import setup
from Cython.Build import cythonize

setup(
        name='integrate app',
        ext_modules=cythonize("integrate.pyx"),
        zip_safe=False,
)
