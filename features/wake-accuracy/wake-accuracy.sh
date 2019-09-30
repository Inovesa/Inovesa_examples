#!/bin/bash

if [ $# -eq 0 ]
then
  inovesa_bin=inovesa
else
  inovesa_bin=$1
fi

mkdir -p results

padding=0002
basename0="results/free-space-pad${padding}"
${inovesa_bin} -c wake-accuracy.cfg -o ${basename0}.png --padding ${padding}

for padding in 0002 0256 0512 1024
do
  basename="results/free-space-pad${padding}"
  ${inovesa_bin} -c wake-accuracy.cfg --padding ${padding} -o ${basename}.h5
  ./plot_wake.py ${basename}.h5
  ${inovesa_bin} -c wake-accuracy.cfg --padding ${padding} -o ${basename}_png.h5 -i ${basename0}.png
  ./plot_wake.py ${basename}_png.h5
done

./wake-accuracy.py
