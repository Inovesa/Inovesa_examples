#! /usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import h5py
import sys

basename = "results/free-space-pad"
paddings = ["1024", "0512", "0004", "0002"]

fig1 = plt.figure(figsize=(4, 3), dpi=80, tight_layout=True)
ax1 = fig1.add_subplot(111)
ax1.set_xlabel(r"Long. Position ($\sigma_q$)")
ax1.set_xlim(-6,6)
ax1.set_ylabel(r"Quantity (arb. units)")
ax1.set_ylim(-1,1)

fname_flt = basename+"1024.h5"
fname_int = basename+"1024_png.h5"
hdf_flt = h5py.File(fname_flt, 'r')
hdf_int = h5py.File(fname_int, 'r')

ax1.plot(hdf_flt["/BunchProfile/axis1"][:],
            hdf_flt["/BunchProfile/data"][0,0,:],
            "k-",
            label="Charge Density (a)")
ax1.plot(hdf_int["/BunchProfile/axis1"][:],
            hdf_int["/BunchProfile/data"][0,0,:],
            "g:",
            label="Charge Density (b)")
ax1.plot(hdf_flt["/BunchProfile/axis1"][:],
            hdf_flt["/WakePotential/data"][0,0,:],
            "r-",
            label="Wake Potential (a)")
ax1.plot(hdf_int["/BunchProfile/axis1"][:],
            hdf_int["/WakePotential/data"][0,0,:],
            "b:",
            label="Wake Potential (b)")

ax1.set_yscale("symlog", linthreshy=1e-9)
ax1.locator_params(axis='y', numticks=7)
ax1.grid()
ax1.legend(loc='center')

fig1.savefig("wake-accuracy-step.eps")


fig2 = plt.figure(figsize=(4, 3), dpi=80, tight_layout=True)
ax2 = fig2.add_subplot(111)
ax2.set_xlabel(r"Long. Position ($\sigma_q$)")
ax2.set_xlim(-6,6)
ax2.set_ylabel(r"Quantity (arb. units)")


ax2.plot(hdf_int["/BunchProfile/axis1"][:],
            hdf_int["/BunchProfile/data"][0,0,:],
            "k:",
            label="Charge Density")

hdf_flt.close()
hdf_int.close()

for padding in paddings:
    fname_flt = basename+padding+".h5"
    fname_int = basename+padding+"_png.h5"
    hdf_flt = h5py.File(fname_flt, 'r')
    hdf_int = h5py.File(fname_int, 'r')
  
    ax2.plot(hdf_flt["/BunchProfile/axis1"][:],
             hdf_int["/WakePotential/data"][0,0,:],
             label="Padding="+padding)

    hdf_flt.close()
    hdf_int.close()

ax2.set_yscale("symlog", linthreshy=1e-6)
ax2.set_ylim(-1,1)
ax2.grid()
ax2.legend(loc='center')

fig2.savefig("wake-accuracy_pad.eps")
