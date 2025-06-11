#Joshua Pan Nov 2023
#Following Allen and Ingram, plot log-log area-weighted CDF of hourly grid-cell precip

#PRECDIR = '/glade/work/jpan/PRECTux'
GRIDDIR = '/glade/p/cesmdata/inputdata/share/scripgrids'
GRIDFN = 'ne120np4_pentagons_100310.nc'

HISTS = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/*.h1i.*.nc'
CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
ALIASES = ['UNSEED', 'CTRL', 'SEED']

cdfs = [] #store cdf of each case
CDFPKL = '/glade/work/jpan/PRECTbmom/prect_cdf_8.pkl'
CALCCDFS = True

a = 6.371e6 #m
MPS2MMHR = 3.6e6
thresh = 1e-8 #m s-1, drizzle threshold

import sys
import os
import numpy as np
import dask.array as da
import uxarray as ux
import xarray as xr
import pickle

import pltsettings
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def main():
   global cdfs
   pltsettings.set1()

   if CALCCDFS:
      for case in CASES:
         cdfs.append(compute_cdfs(HISTS % case))
   
      print('Pickling...')
      with open('/glade/work/jpan/PRECTbmom/prect_cdf_%d.pkl' % int(-np.log10(thresh)), 'wb') as fl:
         pickle.dump(cdfs, fl)
   else:
      with open(CDFPKL, 'rb') as fl:
         cdfs = pickle.load(fl)

   pltcdf_line(cdfs[0], ALIASES[0])
   pltcdf_line(cdfs[1], ALIASES[1], ax=plt.gca())
   pltcdf_line(cdfs[2], ALIASES[2], ax=plt.gca())

   plt.legend()
   plt.savefig('PRECT_cdf_test.png')
   plt.close()

def compute_cdfs(histpath, latbins=None):
   print('Computing CDFs', histpath)
   ds = ux.open_mfdataset(os.path.join(GRIDDIR, GRIDFN), histpath)
   prec1d = ds.PRECT.values.reshape(-1) #shape (time, n_face)
   areas = (ds.uxgrid.face_areas * a**2).data.repeat(ds.time.shape[0])

   print('\tSorting...')
   sorter = np.argsort(prec1d)
   precsort = prec1d[sorter]
   wgts = areas[sorter]
   wherethresh = np.where(precsort > thresh)[0][0]
   precsort = precsort[wherethresh:]
   wgts = wgts[wherethresh:]
   qtiles = np.cumsum(wgts) / sum(wgts)

   return precsort, wgts, qtiles

#takes cdf tuple like compute_cdfs() output
def pltcdf_line(cdftup, label, ax=None, **pltkwargs):
   if ax is None:
      ax = plt.axes()
      ax.set_xscale('log')
      ax.set_yscale('log')
      #ticks = np.logspace(-8, -1, 8)
      #plt.xticks(ticks, labels=np.log10(ticks))
      ticklocs = 10. ** np.arange(-6, 0.1, 2)
      ticklabs = 100 * (1 - ticklocs)
      ax.set_xticks(ticklocs, labels=ticklabs)
      ax.invert_xaxis()
      ax.set_xlabel('Percentile')
      ax.set_ylabel('P rate [mm h$^{-1}$]')

   ln = ax.plot(1 - cdftup[2], cdftup[0] * MPS2MMHR, label=label)
   return ax, ln


   '''
      print('Plotting...')      
      plt.plot(1 - qtiles, precsort * MPS2MMHR, label=ALIASES[ii])
      if case == CASES[0]:
         plt.xscale('log')
         plt.yscale('log')
         #ticks = np.logspace(-8, -1, 8)
         #plt.xticks(ticks, labels=np.log10(ticks))
         ticklocs = 10. ** np.arange(-6, 0.1, 2)
         ticklabs = 100 * (1 - ticklocs)
         plt.xticks(ticklocs, labels=ticklabs)
         plt.gca().invert_xaxis()
         plt.xlabel('Percentile')
         plt.ylabel('P rate [mm h$^{-1}$]')
      if case == CASES[-1]:
         plt.legend(loc='center right')
      plt.savefig('%s_prect_cdf_%d.png' % (case, int(-np.log10(thresh))), bbox_inches='tight')

   print('Computing ratios at percentiles...')
   rattiles = np.logspace(-8, -0.5, 100)
   pp, qq = 2, 1
   rats = []
   for tile in rattiles:
      idx0 = np.where(cdfs[pp][2] >= 1 - tile)[0][0]
      idx1 = np.where(cdfs[qq][2] >= 1 - tile)[0][0]
      p0 = cdfs[pp][0][idx0]
      p1 = cdfs[qq][0][idx1]
      rats.append(p0 / p1)
   ax0 = plt.gca().twinx()
   ax0.plot(rattiles, rats, linestyle='-.', color='black')
   ax0.axhline(1, linestyle='--', color='gray')
   ax0.set_ylabel('%s/%s' % (ALIASES[pp], ALIASES[qq]))
   plt.savefig('prect_cdf_ratios_%d.png' % int(-np.log10(thresh)), bbox_inches='tight')
   plt.close()



   print('%s done.' % sys.argv[0])
   '''

   '''
   print('Plotting all percentiles to find threshold...')
   plt.plot([i[1] for i in oneto100])
   plt.yscale('log')
   plt.savefig('precvol_1to100.png', bbox_inches='tight')
   plt.close()
   '''

if __name__ == '__main__':
   main()
