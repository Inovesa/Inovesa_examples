#! /usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import h5py
import sys
from scipy.interpolate import interp1d

fname = sys.argv[1]
h5file = h5py.File(fname, 'r')

nbl = h5file["/Info/AxisValues_z"].attrs["Second"]*1e12
keV = h5file["/WakePotential/data"].attrs["Volt"]/1e3
steps_per_Ts = h5file["/Info/Parameters"].attrs["steps"]
angle = 2*np.pi/steps_per_Ts
outstep = h5file["/Info/Parameters"].attrs["outstep"]
revolutions_per_step = h5file["/Info/AxisValues_t"][1]*h5file["/Info/AxisValues_t"].attrs["Turn"]/outstep
n_gridpoints = h5file["/Info/Parameters"].attrs["GridSize"]
gridshift_x = h5file["/Info/Parameters"].attrs["PhaseSpaceShiftX"]
edgepoints = [(-(n_gridpoints/2.0)),n_gridpoints/2.0-2*gridshift_x]

rf_flank = np.tan(angle)*np.linspace(edgepoints[0],edgepoints[1],n_gridpoints)*keV
axis_z = np.array(h5file["/Info/AxisValues_z"])*nbl
bunchprof = h5file["/BunchProfile/data"][-1]*h5file["/BunchProfile/data"].attrs["CoulombPerNBL"]*nbl*1e12
wakepot = h5file["/WakePotential/data"][-1]*keV
com = h5file["/BunchPosition/data"][-1]*h5file["/BunchPosition/data"].attrs["Second"]*1e12
h5file.close()

deltaE=(rf_flank-wakepot)/revolutions_per_step
invf = interp1d(deltaE, axis_z)
syncphase = invf(0)




plt.figure(figsize=(3.5,2),tight_layout=True)

plt.plot(axis_z,deltaE,"b")
plt.plot([syncphase,syncphase],[-30,30],"b--")

plt.ylabel("Energy gain (keV)")
plt.ylim(-27.5,27.5)
plt.yticks(np.arange(-20, 30, step=10))

plt.gca().spines['right'].set_color('red')
plt.gca().spines['left'].set_color('blue')
plt.gca().yaxis.label.set_color('blue')
plt.gca().tick_params(axis='y', colors='blue')


plt.grid()


plt.xlabel("Longitudinal Position (ps)")

plt.twinx()

plt.xlabel("Longitudinal position (ps)")
plt.ylabel("Charge density (pC/ps)")
plt.plot(axis_z,bunchprof,"r-")
plt.plot([com,com],[0,15],"r--")

plt.gca().spines['right'].set_color('red')
plt.gca().spines['left'].set_color('blue')
plt.gca().yaxis.label.set_color('red')
plt.gca().tick_params(axis='y', colors='red')

plt.yticks(np.arange(0, 14, step=3))
plt.ylim(0,14)

plt.xticks(np.arange(-6, 8, step=3))
plt.xlim(-8,8)

plt.savefig("plots/f05493_example-syncphase.eps")
plt.savefig("plots/f05493_example-syncphase.pdf")
plt.savefig("plots/f05493_example-syncphase.png")


