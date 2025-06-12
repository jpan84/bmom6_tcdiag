#Joshua Pan Nov 2023
#Following Allen and Ingram, plot log-log area-weighted CDF of hourly grid-cell precip

#PRECDIR = '/glade/work/jpan/PRECTux'
GRIDDIR = '/glade/p/cesmdata/inputdata/share/scripgrids'
GRIDFN = 'ne120np4_pentagons_100310.nc'

HISTS = '/glade/derecho/scratch/jpan/archive/%s/atm/hist_0010_h1i/*.h1i.*.nc'
CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
ALIASES = ['UNSEED', 'CTRL', 'SEED']
COLORS = ['blue', 'orange', 'green']
LATBINS = None
MARKERS = ['.', 'v', '1', 's', '*', '+']

cdfs = [] #store cdf of each case as list(tuple(prate, areawgt, quantile))
cdfslat = [] #store cdf of each case as list(list(tuple(prate, areawgt, quantile, (lowerlat, upperlat))))
#outer list indexed by case, inner list indexed by latitude bin
CDFPKL = '/glade/work/jpan/PRECTbmom/prect_cdf_8.pkl'
CDFPKLLAT = '/glade/work/jpan/PRECTbmom/prect_cdf_8_lat.pkl'
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

LATCDF = bool(int(sys.argv[1]))

def main_global():
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

   ratio_line(cdfs[2], cdfs[0], ALIASES[2], ALIASES[0], plt.gca())
   plt.close()

def main_lat():
   global cdfslat
   pltsettings.set1()

   if CALCCDFS:
      for case in CASES: 
         cdfslat.append(compute_cdfs_lat(HISTS % case))
   
         print('Pickling...')
         with open('/glade/work/jpan/PRECTbmom/prect_cdf_%d_lat.pkl' % int(-np.log10(thresh)), 'wb') as fl:
            pickle.dump(cdfslat, fl)
   else:
      with open(CDFPKLLAT, 'rb') as fl:
         cdfslat = pickle.load(fl)

   ax = plt.axes()
   for casei in range(len(ALIASES)):
      for latj in [6, 7, 8, 10]:
         cdftup = cdfslat[casei][latj]
         pltcdf_line(cdftup, ALIASES[casei] + '_' + cdftup[3][0] + '_' + cdftup[3][1], ax=ax, marker=MARKERS[latj])

   plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=8))
   plt.savefig('PRECT_cdf_lat_test.png')
   plt.close()

def compute_cdfs(histpath):
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

def compute_cdfs_lat(histpath, latbins=np.arange(-90, 90.1, 15)):
   print('Computing lat-binned CDFs', histpath)
   ds = ux.open_mfdataset(os.path.join(GRIDDIR, GRIDFN), histpath)

   cdftups = []#precsort, wgts, qtiles = None, None, None
   for ii in range(latbins.size - 1):
      llat, ulat = latbins[ii], latbins[ii + 1]
      print('\t lats %.1f %.1f' % (llat, ulat))
      idx = np.where((ds.lat >= llat) & (ds.lat < ulat))[0]

      prec1d = ds['PRECT'].isel(n_face=idx).data.reshape(-1) #shape (time, n_face)
      areas = (ds['area'].isel(n_face=idx) * a**2).data.repeat(ds.time.size, axis=0)

      print('\tSorting...')
      sorter = np.argsort(prec1d)
      precsort = prec1d[sorter]
      wgts = areas[sorter]
      wherethresh = np.where(precsort > thresh)[0][0]
      precsort = precsort[wherethresh:]
      wgts = wgts[wherethresh:]
      qtiles = np.cumsum(wgts) / sum(wgts)

      cdftups.append((precsort, wgts, qtiles, (llat, ulat)))

   return cdftups

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

   ln = ax.plot(1 - cdftup[2], cdftup[0] * MPS2MMHR, label=label, **pltkwargs)
   return ax, ln

def ratio_line(cdftup1, cdftup2, alias1, alias2, ax):
   print('Computing ratios for %s / %s...' % (alias1, alias2))
   rattiles = np.logspace(-8, -0.5, 100)
   rats = []
   for tile in rattiles:
      idx1 = np.where(cdftup1[2] >= 1 - tile)[0][0]
      idx2 = np.where(cdftup2[2] >= 1 - tile)[0][0]
      p1 = cdftup1[0][idx1]
      p2 = cdftup2[0][idx2]
      rats.append(p1 / p2)
   ax1 = ax.twinx()
   ax1.plot(rattiles, rats, linestyle='-.', color='black')
   ax1.axhline(1, linestyle='--', color='gray')
   ax1.set_ylabel('%s/%s' % (alias1, alias2))
   plt.savefig('prect_cdf_ratios_%d.png' % int(-np.log10(thresh)), bbox_inches='tight')
   #plt.close()

if __name__ == '__main__':
   if LATCDF:
      main_lat()
   else:
      main_global()
