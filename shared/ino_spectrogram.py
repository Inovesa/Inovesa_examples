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
  parser = argparse.ArgumentParser(description='Generating CSR Spectrogramms from Inovesa result files')
  parser.add_argument('directory', type=str, help='relativ path of the Inovesa hdf5 files')
  parser.add_argument('--ending', type=str, default='b.h5', help='only files ending with the specified string are used')
  parser.add_argument('--saveplot', type=str, nargs='?', default=False, const='replace by args.directory', help='save plot under filename (including path). Without additional argument, plot is save with standard name in path of hdf5 files. ')
  parser.add_argument('--showplot', action='store_true', help='shows plot')
  parser.add_argument('--datalen', type=int, default=0, help='spezifies how many simulations steps (starting from the end of the file) will be used')
  parser.add_argument('--currentmin', type=float,default=None, help='minimum of bunch current axis')
  parser.add_argument('--currentmax', type=float,default=None, help='maximum of bunch current axis')
  parser.add_argument('--freqmin', type=float,default=None, help='minimum of frequency axis')
  parser.add_argument('--freqmax', type=float,default=None, help='maximum of frequency axis')
  parser.add_argument('--colmin', type=float,default=None, help='minimum of intensity (color) axis')
  parser.add_argument('--colmax', type=float,default=None, help='maximum of intensity (color) axis')
  parser.add_argument('--title', action='store_true', help='print title in plot (derived from directory)')
  parser.add_argument('--scaling', type=float,default=1, help='scaling factor for conversion from isomagnetic to real ring (if needed)')
  parser.add_argument('--xlog', action='store_true', help='set x axis logarithmic')
  parser.add_argument('--ylog', action='store_true', help='set y axis logarithmic')

  args = parser.parse_args()
  args.saveplot = "%s/simulated-spectogram.png" %(args.directory) if args.saveplot=='replace by args.directory' else args.saveplot

  if not args.showplot and not args.saveplot:
    print("Nothing to do?!?")

  import matplotlib

  if not args.showplot:
      matplotlib.use('Agg')
  if args.saveplot:
    if args.saveplot.split('.')[-1] == 'eps':
      matplotlib.use('PS')

  import matplotlib.pyplot as plt
  from matplotlib.ticker import FuncFormatter

  dataname='/CSR/Intensity/data'
  timename='/Info/AxisValues_t'

  fnames = []
  
  for fname in os.listdir(args.directory):
    if fname[-len(args.ending):] == args.ending:
      fnames.append(fname)
  
  if len(fnames) == 0:
    print("No such file: '"+ args.directory+"/*"+args.ending+"'")
    exit()

  title = filter(None,args.directory.split('/'))[-1] if args.title else ""

  data = []
  currents = []
  deltat = 1
  
  datalen = args.datalen

  for fname in fnames:
    try:
      h5file1 = h5py.File(args.directory+'/'+fname, 'r')
      tmp=deltat
      try: # Inovesa v0.14
          deltat = h5file1['/Info/AxisValues_t'].attrs['Second']*h5file1['/Info/AxisValues_t'][1]
      except: # older versions
          deltat = h5file1['/Info/AxisValues_t'].attrs['Factor4Seconds']*h5file1['/Info/AxisValues_t'][1]
      assert tmp==deltat or tmp==1, 'not same deltat in all files (at file %s)' %fname
      if datalen == 0:
      	  datalen = int(len(h5file1[dataname][...])/2-1)
      assert len(h5file1[dataname][...]) > 2*datalen, 'to few data points (in file %s)' %fname
      data.append(np.abs(np.fft.rfft(h5file1[dataname][-2*datalen:]))[1:])
      currents.append(h5file1['/Info/Parameters'].attrs['BunchCurrent'])
      h5file1.close()
    except IOError as e:
      print(e, 'skipping %s and continuing with next file' %fname)

  freqs = np.linspace(0,args.scaling/(2*deltat),datalen)/1000.

  data = np.array(data)
  currents = np.array(currents)
  currents = currents*1000
  data = (data.T*np.power(currents,2)).T
  data=data[currents.argsort(),:]
  currents.sort()
  currents = currents*args.scaling

  if args.saveplot:
    if args.saveplot.split('.')[-1]=='npz':
      npzfname = os.path.dirname(args.saveplot)+os.path.basename(args.saveplot)
      np.savez(npzfname,freqs=freqs,currents=currents,data=data)
      print("Saving data to " + args.saveplot)
      exit()

  font = {'size' : 19 }
  matplotlib.rc('font', **font)

  plt.figure(figsize=(10,5),tight_layout=True)
  plt.xlabel("Frequency / %.2f kHz" %args.scaling if args.scaling !=1 else "Frequency (kHz)")
  plt.ylabel("Bunch Current / %.2f mA" %args.scaling if args.scaling !=1 else "Bunch Current (mA)")
  spectrogram = plt.pcolormesh(freqs,currents,data,shading='flat',norm=matplotlib.colors.LogNorm(vmin=args.colmin,vmax=args.colmax))
  if args.ylog:
    plt.yscale('log', subsy=np.arange(1, 10, 1)[1:])
  plt.gca().yaxis.set_minor_formatter(FuncFormatter(minor_format))
  plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, i=None: x))
  plt.tick_params('x',which='major', direction='out', labelleft='on', width=2,length=5)
  plt.tick_params('y',which='minor', direction='inout', labelleft='on', width=1,length=5)
  if args.xlog:
    plt.xscale('log', subsy=np.arange(1, 10, 1)[1:])
  plt.gca().yaxis.set_ticks_position('left')
  plt.gca().xaxis.set_ticks_position('bottom') 
  plt.xlim((args.freqmin if args.freqmin else 0.07 if args.xlog else freqs[0], args.freqmax if args.freqmax else freqs[-1]))
  plt.ylim((args.currentmin if args.currentmin else currents[0],args.currentmax if args.currentmax else currents[-1]))
  if args.ylog:
    plt.minorticks_on()
  plt.colorbar(spectrogram).set_label("Spectral Intensity (arb. units)")
  spectrogram.set_cmap('inferno')
  plt.title(title)

  if args.saveplot:
    print("Saving plot to " + args.saveplot)
    plt.savefig(args.saveplot,dpi=200)

  if args.showplot:
    plt.show()
  else:
    plt.close()


if __name__ == '__main__':

  main()
