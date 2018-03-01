#! /bin/bash

# Change, e.g. if you have not installed Inovesa.
#binary="/home/user/inovesa/inovesa"
#binary="inovesa"
binary="/home/schoenfeldt/Projects/build/inovesa/release/inovesa"

# We create a subdirectory named after the version of Inovesa.
inoversion=`$binary --version`
inoversion=${inoversion%%,*}
inoversion=${inoversion/ /_}
dir="./${inoversion}"

ending="_a.h5"

mkdir -p plots

for f in ${dir}/*_b.h5
do
  if [ -e "$f" ]
  then
    ending="_b.h5"
  fi
  
  # we just used the loop as a tool to find whether there is at least one file
  break
done

../../shared/ino_spectrogram.py --saveplot plots/f05135_spectrogram.png ${dir} --ending ${ending} --scaling 0.31621
