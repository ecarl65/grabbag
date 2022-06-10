#!/bin/bash

debug="-g"
if [[ $# -ge 1 && $1 -ne "0" ]]; then
        debug+=" -DDEBUG"
fi

downsamp="8"
if [[ $# -ge 2 ]]; then
        downsamp=$2
fi

filtord=8
if [[ $# -ge 3 ]]; then
        filtord=$3
fi

make clean
make CFLAGS="${debug}" fudft
./fudft -d${downsamp} -f${filtord}
./plot_fudft.py ${downsamp}
