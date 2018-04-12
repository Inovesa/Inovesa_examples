#! /usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import LogLocator
import numpy as np
import os
import h5py
import sys
import scipy.signal
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

def minor_format(x, i=None):
    return x

def fwhm(data):
    fwhm = np.zeros(np.shape(data)[0])
    for i,profile in enumerate(data):
        maxval = max(profile)
        maxpos = np.where(profile==maxval)[0]
        above = np.size((np.where(profile>maxval/2)))
        fwhm[i] = above
    return fwhm

fnames = []

ending = "_a.h5"
for fname in os.listdir(sys.argv[1]):
  if len(fname) == len("1234_a.h5") and fname[-len(ending):] == ending:
    fnames.append("./v1.0.0/"+fname)

datapoints = 2000
dataset = []

zerobunchlength = 1

for fname in fnames:
  try:
    h5file1 = h5py.File(fname, 'r')
    
    current = (h5file1['/Info/Parameters'].attrs['BunchCurrent'])
    
    nbl = h5file1["/Info/AxisValues_z"].attrs["Second"]*1e12
    keV = h5file1["/WakePotential/data"].attrs["Volt"]/1e3
    angle = 2*np.pi/h5file1["/Info/Parameters"].attrs["steps"]
    n_gridpoints = h5file1["/Info/Parameters"].attrs["GridSize"]
    gridshift_x = h5file1["/Info/Parameters"].attrs["PhaseSpaceShiftX"]
    edgepoints = [(-(n_gridpoints/2.0)),n_gridpoints/2.0-2*gridshift_x]

    rf_flank = np.tan(angle)*np.linspace(edgepoints[0],edgepoints[1],n_gridpoints)*keV
    axis_z = np.array(h5file1["/Info/AxisValues_z"])*nbl
    wakepot = h5file1["/WakePotential/data"][:datapoints,]*keV
    bunchposition = h5file1["/BunchPosition/data"][:datapoints,]*h5file1["/BunchPosition/data"].attrs["Second"]*1e12
    
    csrloss = h5file1["/CSR/Intensity/data"][:datapoints,]*h5file1["/CSR/Intensity/data"].attrs["Watt"]
    
    syncphase = []
    deltaE=rf_flank-wakepot
    for i in range(len(bunchposition)):
        invf = interp1d(deltaE[i,:], axis_z)
        syncphase.append(invf(0))
    
    syncphase = np.array(syncphase)
    
    dataset.append([current,csrloss,bunchposition,syncphase])
    h5file1.close()
  except:
    print(fname)


dataset.sort()

currents = [ 0 ]

csrloss_mean = [ 0 ]
csrloss_min = [ 0 ]
csrloss_max = [ 0 ]
csrloss_sdev = [ 0 ]
csrloss_sdev_lower = [ 0 ]
csrloss_sdev_upper = [ 0 ]

bunchposition_mean = [ 0 ]
bunchposition_min = [ 0 ]
bunchposition_max = [ 0 ]
bunchposition_sdev = [ 0 ]
bunchposition_sdev_lower = [ 0 ]
bunchposition_sdev_upper = [ 0 ]

syncphase_mean = [ 0 ]
syncphase_min = [ 0 ]
syncphase_max = [ 0 ]
syncphase_sdev = [ 0 ]
syncphase_sdev_lower = [ 0 ]
syncphase_sdev_upper = [ 0 ]

blDelta = [ 0 ]

for current,csrloss,bunchposition,syncphases in dataset:
    currents.append(current)
    
    csrloss = csrloss/current/1e3
    csrloss_mean.append(np.mean(csrloss))
    csrloss_sdev.append(np.std(csrloss))
    mask = csrloss < csrloss_mean[-1]
    csrloss_sdev_upper.append(np.ma.masked_where(mask, csrloss).std())
    csrloss_sdev_lower.append(np.ma.masked_where(~mask, csrloss).std())
    csrloss_min.append(np.min(csrloss))
    csrloss_max.append(np.max(csrloss))
    
    bunchposition_mean.append(np.mean(bunchposition))
    bunchposition_sdev.append(np.std(bunchposition))
    mask = bunchposition < bunchposition_mean[-1]
    bunchposition_sdev_upper.append(np.ma.masked_where(mask, bunchposition).std())
    bunchposition_sdev_lower.append(np.ma.masked_where(~mask, bunchposition).std())
    bunchposition_min.append(np.min(bunchposition))
    bunchposition_max.append(np.max(bunchposition))
    
    syncphase_mean.append(np.mean(syncphases))
    syncphase_sdev.append(np.std(syncphases))
    mask = syncphases < syncphase_mean[-1]
    syncphase_sdev_upper.append(np.ma.masked_where(mask, syncphases).std())
    syncphase_sdev_lower.append(np.ma.masked_where(~mask, syncphases).std())
    syncphase_min.append(np.min(syncphases))
    syncphase_max.append(np.max(syncphases))

currents = np.array(currents)*1e3

csrloss_mean = np.array(csrloss_mean)
csrloss_min = np.array(csrloss_min)
csrloss_max = np.array(csrloss_max)
csrloss_sdev = np.array(csrloss_sdev)

bunchposition_mean = np.array(bunchposition_mean)
bunchposition_min = np.array(bunchposition_min)
bunchposition_max = np.array(bunchposition_max)
bunchposition_sdev = np.array(bunchposition_sdev)

syncphase_mean = np.array(syncphase_mean)
syncphase_min = np.array(syncphase_min)
syncphase_max = np.array(syncphase_max)
syncphase_sdev = np.array(syncphase_sdev)


#plt.rc('text', usetex=True)
#font = {'family' : 'normal',  'weight' : 'bold', 'size' : 16 }
#font = {'size' : 19 }
#matplotlib.rc('font', **font)


plt.figure(figsize=(3.5,2),tight_layout=True)

plt.xlabel("Bunch Current (mA)")
plt.ylabel(r"$\phi(\Delta E = 0)$ (ps)")
plt.gca().spines['right'].set_color('red')
plt.gca().spines['left'].set_color('blue')
plt.gca().yaxis.label.set_color('blue')
plt.gca().tick_params(axis='y', colors='blue')

#plt.fill_between(currents,syncphase_min,syncphase_max,color="#FFDDDD",label="",hatch='/')
plt.fill_between(currents,syncphase_mean-syncphase_sdev_lower,syncphase_mean+syncphase_sdev_upper, color="#AAAAFF")
handle0, = plt.plot(currents,syncphase_mean,"b-",label=r"$\phi(\Delta E = 0)$")

plt.xscale("log")
plt.ylim(-2,3)
plt.yticks(np.arange(0, 4, step=1))

ax1 = plt.gca()

plt.twinx()


plt.xlabel("Bunch Current (mA)")
plt.ylabel(r"$\langle \phi \rangle$ (ps)")
plt.gca().spines['right'].set_color('red')
plt.gca().spines['left'].set_color('blue')
plt.gca().yaxis.label.set_color('red')
plt.gca().tick_params(axis='y', colors='red')

#plt.fill_between(currents,bunchposition_min,bunchposition_max,facecolor="none",edgecolor="#DDDDFF",label="",hatch='\\')
plt.fill_between(currents,bunchposition_mean-bunchposition_sdev_lower,bunchposition_mean+bunchposition_sdev_lower,facecolor="r",label="",alpha=0.5, linewidth=0.0)
handle1, = plt.plot(currents,bunchposition_mean,"r--",label=r"$\langle \phi \rangle$")

plt.xlim(0.001,1)
plt.ylim(0,1)
plt.yticks(np.arange(0, 1, step=0.2))
plt.grid()



plt.savefig("plots/f05493_syncphase-shift.eps")
plt.savefig("plots/f05493_syncphase-shift.pdf")
plt.savefig("plots/f05493_syncphase-shift.png")

