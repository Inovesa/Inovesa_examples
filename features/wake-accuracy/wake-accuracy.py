#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import h5py
import sys

fname = sys.argv[1]

hdf_f = h5py.File(fname, 'r')

axis_x = hdf_f["/BunchProfile/axis1"][:]
bunchprofile = hdf_f["/BunchProfile/data"][0,0,:]
wakepotential = hdf_f["/WakePotential/data"][0,0,:]
inovesa_version = hdf_f['/Info/Inovesa_v'][:]

if (inovesa_version[:2] != [1,1]).any():
    print("Only tested for Inovesa 1.1.")

hdf_f.close()

plt.figure(figsize=(6, 4), dpi=80, tight_layout=True)
ax1 = plt.subplot(111)
ax1.set_xlabel(r"Long. Position ($\sigma_q$)")
ax1.set_xlim(-6,6)
ax1.set_ylabel(r"Quantity (arb. units)")
ax1.set_ylim(-1,1)

ax1.plot(axis_x, bunchprofile, label="Charge Density")
ax1.plot(axis_x, wakepotential, label="Wake Potential")

ax1.set_yscale("symlog", linthreshy=1e-9)
ax1.grid()
ax1.legend()

plt.savefig(fname+".eps")

