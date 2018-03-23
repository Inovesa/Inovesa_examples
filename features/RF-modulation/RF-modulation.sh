#! /bin/bash

# Change, e.g. if you have not installed Inovesa.
#binary="/home/user/inovesa/inovesa"
binary="inovesa"
binary="/home/patrik/Projekte/build/inovesa/release/inovesa"

# We create a subdirectory named after the version of Inovesa.
inoversion=`$binary --version`
inoversion=${inoversion%%,*}

inoversion_major=${inoversion%%.*}
inoversion_major=${inoversion_major#v}
inoversion_minor=${inoversion%% *}
inoversion_minor=${inoversion_minor#*.}
inoversion_minor=${inoversion_minor%.*}

if [[ $inoversion_major -le 1 && $inoversion_minor -le 0 ]]
then
  echo "Found Inovesa ${inoversion}."
  echo "But need at least v1.1."
  exit
fi


inoversion=${inoversion/ /_}
dir="./${inoversion}"
mkdir -p $dir

example="RF-modulation_slow-phase"

$binary -c ${example}.cfg -o ${dir}/$example.h5

echo "Plotting the results..."

python3 ../../shared/plot-hdf5.py mq ${dir}/$example.h5



