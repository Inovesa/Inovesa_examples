#! /usr/bin/env python3
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
import io
import unicodecsv

meas_file = "f06112_meas.csv"

timestamp = []
meas_current34 = []
meas_current35 = []
meas_timing34 = []
meas_timing35 = []
acc_timing34 = []
acc_timing35 = []

with io.open(meas_file, encoding='utf8') as f:
	reader = unicodecsv.reader(f,delimiter=',')
	for entry in reader:
		if entry[0][0] != 'U':
			timestamp.append(int(entry[1]))
			meas_current34.append(float(entry[5]))
			meas_current35.append(float(entry[6]))
			meas_timing34.append(float(entry[7]))
			acc_timing34.append(float(entry[8]))
			meas_timing35.append(float(entry[9]))
			acc_timing35.append(float(entry[10]))


timestamp = np.array(timestamp)
meas_timing34 = np.array(meas_timing34)
meas_timing35 = np.array(meas_timing35)
acc_timing34 = np.array(acc_timing34)
acc_timing35 = np.array(acc_timing35)
delta_errorbar = np.sqrt(acc_timing34*acc_timing34+acc_timing35*acc_timing35)


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
plt.ylabel(r"$\langle \phi \rangle$ (ps)")

#plt.fill_between(currents,bunchposition_min,bunchposition_max,facecolor="none",edgecolor="#DDDDFF",label="",hatch='\\')
plt.fill_between(currents,bunchposition_mean-bunchposition_sdev_lower,bunchposition_mean+bunchposition_sdev_lower,facecolor="m",label="",alpha=0.5, linewidth=0.0)
handle1, = plt.plot(currents,bunchposition_mean,"m-",label=r"$\langle \phi \rangle$")

plt.xlim(0.001,1.5)
plt.xscale("log")
plt.yticks(np.arange(0, 0.6, step=0.1))
plt.ylim(0,0.5)
plt.grid()



plt.savefig("plots/f06112_syncphase.eps")
plt.savefig("plots/f06112_syncphase.pdf")
plt.savefig("plots/f06112_syncphase.png")

plt.close()


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
plt.yticks(np.arange(0, 1.6, step=0.5))
plt.ylim(-1,1.5)

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

plt.xlim(0.001,1.5)
plt.yticks(np.arange(0, 0.6, step=0.1))
plt.ylim(0,0.5)
plt.grid()



plt.savefig("plots/f06112_syncphase-shift.eps")
plt.savefig("plots/f06112_syncphase-shift.pdf")
plt.savefig("plots/f06112_syncphase-shift.png")

plt.close()



I1 = interp1d(timestamp, meas_current34)
I2 = interp1d(timestamp, meas_current35)
phi = interp1d(currents, bunchposition_mean)

plt.figure(figsize=(3.5,3),tight_layout=True)

gs = matplotlib.gridspec.GridSpec(2, 1)

ax1b = plt.subplot(gs[0])
ax1a = ax1b.twinx()
ax2 = plt.subplot(gs[1])

ax1a.set_ylabel(r"$I$ (mA)")
ax1a.spines['right'].set_color('blue')
ax1a.spines['left'].set_color('red')
ax1a.yaxis.label.set_color('blue')
ax1a.tick_params(axis='y', colors='blue')

ax1a.plot(timestamp/3600.0,I1(timestamp), "b-", label="$I_1$")
ax1a.plot(timestamp/3600.0,I2(timestamp), "b--", label="$I_2$")

ax1a.set_ylim(0,1.75)



ax1b.set_ylabel(r"$\langle \phi \rangle$ (ps)")
ax1b.spines['right'].set_color('blue')
ax1b.spines['left'].set_color('red')
ax1b.yaxis.label.set_color('red')
ax1b.tick_params(axis='y', colors='red')
ax1b.set_yticks(np.arange(0, 0.4, step=0.1))

plt.setp( ax1b.get_xticklabels(), visible=False)

ax1b.plot(timestamp/3600.0,phi(I1(timestamp)),"r-",label=r"$\langle\phi\rangle_1$")
ax1b.plot(timestamp/3600.0,phi(I2(timestamp)),"r--",label=r"$\langle\phi\rangle_2$")

ax1b.set_xlim(0,8)
ax1b.set_ylim(0,0.35)

ax1a.set_yticks(np.arange(0, 2, step=0.5))
ax1b.grid()

meas_delta = meas_timing34-meas_timing35

ax2.errorbar(timestamp/3600.0,meas_delta,yerr=delta_errorbar,fmt='.',color="#FFAAFF",label=r"$\Delta\langle\phi\rangle$",zorder=0)
ax2.plot(timestamp/3600.0,meas_delta,'m.',label=r"$\Delta\langle\phi\rangle$",markersize=2,alpha=0.5,zorder=1)
N = 50
average_delta = np.convolve(meas_delta, np.ones((N,))/N, mode='valid')
ax2.plot(timestamp[N/2-1:-N/2]/3600.0,average_delta,"m-",label=r"$\Delta\langle\phi\rangle$",zorder=2)
ax2.plot(timestamp/3600.0,phi(I2(timestamp))-phi(I1(timestamp)),"r-",label=r"$\Delta\langle\phi\rangle$",zorder=3)

ax2.set_xlim(0,8)
ax2.set_yticks(np.arange(0, 4, step=1))
ax2.set_ylim(-0.5,4)

ax2.set_xlabel("Time (h)")
ax2.set_ylabel(r"$\Delta\langle \phi \rangle$ (ps)")

ax2.grid()

plt.savefig("plots/f06112_syncphase-compare.pdf")
plt.close()

plt.figure(figsize=(3.5,2),tight_layout=True)

plt.plot(I2(timestamp)-I1(timestamp),phi(I2(timestamp))-phi(I1(timestamp)),"k-")

plt.xlabel("Bunch Current Difference (mA)")
plt.ylabel(r"$\Delta\langle \phi \rangle$ (ps)")

plt.savefig("plots/f06112_delta-syncphase.eps")
plt.savefig("plots/f06112_delta-syncphase.pdf")
plt.savefig("plots/f06112_delta-syncphase.png")



exit()


