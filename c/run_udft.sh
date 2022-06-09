#!/bin/sh

make clean
make CFLAGS='-g' udft
./udft
./plot_udft.py
