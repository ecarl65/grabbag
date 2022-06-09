#!/bin/sh

make clean
make CFLAGS='-g -DDEBUG' udft
./udft
./plot_udft.py
