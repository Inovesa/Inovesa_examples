#! /usr/bin/env python

import numpy as np
import os
import h5py
import sys
import scipy.signal
import argparse

def minor_format(x, i=None):
    return x

def main():
  import matplotlib
  import matplotlib.pyplot as plt
  from matplotlib.ticker import FuncFormatter

  dataname='/ACST_Spectrogram'
  currentname='/fit_w_current2_L'
  
  h5file1 = h5py.File("ACST_spectrogram.h5", 'r')
  data1 = (h5file1[dataname][:10000,:]).T
  currents1 = (h5file1[currentname][1:])*1e6
  freqs1 = np.linspace(0,15.2588*10,1e4)
  h5file1.close()
  
  file2 = np.load("f06214_spectrogram-noise_000_000.npz")
  data2 = file2["data"]
  currents2 = (file2["currents"])*1e3
  freqs2 = file2["freqs"]
  
  file3 = np.load("f06214_spectrogram-noise_001_001.npz")
  data3 = file3["data"]
  currents3 = (file3["currents"])*1e3
  freqs3 = file3["freqs"]

  fig1 = plt.figure(figsize=(8,2.2))
  
  spec = matplotlib.gridspec.GridSpec(ncols=3, nrows=1,wspace=0.05,bottom=0.2,left=0.08,right=0.99)
  f1_ax1 = fig1.add_subplot(spec[2])
  f1_ax1.set_title("Measurement")
  
  f1_ax2 = fig1.add_subplot(spec[0])
  f1_ax2.set_title("Simulation (no noise)")
  
  f1_ax3 = fig1.add_subplot(spec[1])
  f1_ax3.set_title("Simulation (with noise)")
  
  plt.sca(f1_ax1)
  plt.tick_params("y",labelleft="off")
  
  plt.sca(f1_ax3)
  plt.tick_params("y",labelleft="off")
  
  
  f1_ax2.set_ylabel("Bunch Current ($\mu$A)")
  
  f1_ax1.set_xlabel("Frequency (kHz)")
  spectrogram1 = f1_ax1.pcolormesh(freqs1,currents1,data1,shading='flat',norm=matplotlib.colors.Normalize(-120,-40))
  f1_ax1.set_xlim(0,90)
  f1_ax1.set_ylim(175,300)
  spectrogram1.set_cmap('inferno')
  
  f1_ax2.set_xlabel("Frequency (kHz)")
  spectrogram2 = f1_ax2.pcolormesh(freqs2,currents2,data2,shading='flat',norm=matplotlib.colors.LogNorm(vmin=5e-4,vmax=1e2))
  f1_ax2.set_xlim(0,90)
  f1_ax2.set_ylim(175,300)
  spectrogram2.set_cmap('inferno')
  
  f1_ax3.set_xlabel("Frequency (kHz)")
  spectrogram3 = f1_ax3.pcolormesh(freqs3,currents3,data3,shading='flat',norm=matplotlib.colors.LogNorm(vmin=5e-4,vmax=1e2))
  f1_ax3.set_xlim(0,90)
  f1_ax3.set_ylim(175,300)
  spectrogram3.set_cmap('inferno')
  
  
  plt.savefig("f06214_spectrogram-compare.png",dpi=200)


if __name__ == '__main__':
  main()
