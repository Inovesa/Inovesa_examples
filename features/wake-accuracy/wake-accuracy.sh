#!/bin/bash

if [ $# -eq 0 ]
then
  inovesa_bin=inovesa
else
  inovesa_bin=$1
fi

mkdir results

padding=0000
basename0="results/free-space-pad${padding}"
${inovesa_bin} -c wake-accuracy.cfg -o ${basename0}.png --padding ${padding}

for padding in 0000 0016 0032 0064 0128 0256 0512 1024
do
  basename="results/free-space-pad${padding}"
  ${inovesa_bin} -c wake-accuracy.cfg --padding ${padding} -o ${basename}.h5
  ${inovesa_bin} -c wake-accuracy.cfg --padding ${padding} -o ${basename}_png.h5 -i ${basename0}.png
done

