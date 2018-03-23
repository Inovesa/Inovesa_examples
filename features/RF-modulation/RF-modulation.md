## About ##

The scripts contained in this directory give simple examples
what can be done using the RF modulation freature.
When you have **Inovesa v1.1 or newer** installed, you can directly run

    ./RF-modulation.sh

If Inovesa is not installed you have to define the path inside the script.
Notice: Older versions of Inovesa do not support the RF modulation feature.


### RF-modulation_slow-phase ###

This configuration models a slow modulation of the RF phase.
Its frequency is set to one tenth of the synchrtron frequency,
so the bunch center of mass can follow the modulation.
The set value for RFPhaseModAmplitude is 0.5°,
which corresponds to 3.1 ps.
This is roughly twice the value of the bunch length (1.5 ps).

