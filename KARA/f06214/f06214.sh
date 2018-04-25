#! /bin/bash

# Change, e.g. if you have not installed Inovesa.
#binary="/home/user/inovesa/inovesa"
binary="inovesa"

# We create a subdirectory named after the version of Inovesa.
inoversion=`$binary --version`
inoversion=${inoversion%%,*}
echo "Using Inovesa $inoversion"

inoversion=${inoversion/ /_}
dir="./${inoversion}"
mkdir -p ${dir}

# configuration to use
cfg="./f06214.cfg"

${binary} -I ${lasti}e-6 --config $config -T 1000 -o ${dir}/${lasti}_a1.h5
for curri in {0350..0175..1}
do
  if [ -f ${dir}/noise_000_000/${curri}_a1.h5 ]; then
    echo "${dir}/noise_000_000/${curri}_a1.h5 exists." 
  else
    ${binary} -I ${curri}e-6 --verbose -i ${dir}/noise_000_000/${lasti}_a1.h5 -c $config -T 500 -o ${dir}/${curri}_a1.h5
  fi
  if [ -f ${dir}/noise_000_000/${curri}_b1.h5 ]; then
    echo "${dir}/noise_000_000/${curri}_b1.h5 exists." 
  else
    ${binary} -I ${curri}e-6 --verbose -i ${dir}/noise_000_000/${curri}_a1.h5 -c $config -T 200 --outstep 10 -o ${dir}/noise_000_000/${curri}_b1.h5 --RFAmplitudeSpread 0 --RFPhaseSpread 0
  fi
  if [ -f ${dir}/noise_001_001/${curri}_b1.h5 ]; then
    echo "${dir}/noise_001_001/${curri}_b1.h5 exists." 
  else
    ${binary} -I ${curri}e-6 --verbose -i ${dir}/noise_000_000/${curri}_a1.h5 -c $config -T 200 --outstep 10 -o ${dir}/noise_001_001/${curri}_b1.h5 --RFAmplitudeSpread 0.01 --RFPhaseSpread 0.01
  fi
  lasti=$curri
done


./ino_spectrogram.py --saveplot f06214_spectrogram-noise_001_001.npz ${dir}/noise_001_001/ --ending _b1.h5 --currentmin 0.175 --currentmax 0.35  --freqmin 0 --freqmax 100
./ino_spectrogram.py --saveplot f06214_spectrogram-noise_000_000.npz ${dir}/noise_000_000/ --ending _b1.h5 --currentmin 0.175 --currentmax 0.35  --freqmin 0 --freqmax 100

python3 f06214_spectrogram-compare.py

