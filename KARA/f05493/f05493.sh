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
mkdir -p $dir

# configuration to use
cfg="./f05493.cfg"

# Start with highest current (here in uA, see "e-6" below)
for curri in {1000..0110..0010}; do
currents=`echo ${currents} ${curri}`
done
for curri in {0100..0001..0001}; do
currents=`echo ${currents} ${curri}`
done

# We use some initial extra time in the beginning
# to find physical conditions.
# The result will be overwritten in the first iteration of the loop.
lasti=${currents%% *}
$binary -c $cfg -I ${lasti}e-6 -T 500 -o ${dir}/${lasti}_a.h5

for curri in $currents; do
  # Commands are echoed to improve understanding what happens.
  echo ""
  echo "$binary -c $cfg -I ${curri}e-6 -i ${dir}/${lasti}_a.h5 -T 200 -o ${dir}/${curri}_a.h5"
  
  # Simulation using the "adiabatic" method:
  # Starting distribution is always read in from the last run.
  $binary -c $cfg -I ${curri}e-6 -i ${dir}/${lasti}_a.h5 -T 200 -o ${dir}/${curri}_a.h5
  
  # Save current for next loop iteration.
  lasti=$curri
done

exit

