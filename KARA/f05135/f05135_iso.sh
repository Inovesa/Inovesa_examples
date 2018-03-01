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
cfg="./f05135_iso.cfg"

# Start with highest current (here in uA, see "e-6" below)
lasti=2180

# We use some initial extra time in the beginning
# to find physical conditions.
# The result will be overwritten in the first iteration of the loop.
$binary -c $cfg -I ${lasti}e-6 -T 50 -o ${dir}/${lasti}_a.h5


for curri in {2180..0480..0010}; do
  # Commands are echoed to improve understanding what happens.
  echo ""
  echo "$binary -c $cfg -I ${curri}e-6 -i ${dir}/${lasti}_a.h5 -T 200 -o ${dir}/${curri}_a.h5"
  
  # Simulation using the "adiabatic" method:
  # Starting distribution is always read in from the last run.
  $binary -c $cfg -I ${curri}e-6 -i ${dir}/${lasti}_a.h5 -T 200 -o ${dir}/${curri}_a.h5
  
  # If you are not happy with the simple "adiabatic" result
  # you can do a second simulation run that continues the previous data.
  #echo "$binary -c $cfg -I $curri -i ${dir}/${currt}_a.h5 -T 200 -o ${dir}/${curri}_b.h5"
  #$binary -c $cfg -I $curri -i ${dir}/${currt}_a.h5 -T 200 -o ${dir}/${curri}_b.h5
  
  # Save current for next loop iteration.
  lasti=$curri
done

exit

