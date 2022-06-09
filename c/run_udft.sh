#!/bin/bash

downsamp="8"
if [[ $# -ge 1 ]]; then
        downsamp=$1
fi

filtord=8
if [[ $# -ge 2 ]]; then
        filtord=$2
fi

make clean
make CFLAGS='-g' udft
./udft -d ${downsamp} -f ${filtord}
./plot_udft.py ${downsamp}
