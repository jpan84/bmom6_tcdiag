#Joshua Pan Mar 2025
#Uses parquet output of traj_stats_bmom_SE.py

import os
import sys
import pickle
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats

MS2KT = 1.9438
a_e = 6.371e6 #m
LATLABS = np.arange(-90, 90.1, 30)
FIXYLIM = (0, 70)
DLAT = 0.5

#FILI = sys.argv[1]
FILI = ['250415_unseed.parquet', '250417_ctrl.parquet', '250416_seed1x1.parquet']
labels = ['unseed', 'ctrl', 'seed']
DOUT = '250519_density'

#FILI = ['250415_unseed_JJASOND.parquet', '250417_ctrl_JJASOND.parquet']
#labels = ['unseed', 'ctrl']
#DOUT = '250423_density_JJASOND'

def main():
   #extidx = FILI.rfind('.')
   #DOUT = FILI[:extidx] + '_density'
   if not os.path.exists(DOUT):
      os.makedirs(DOUT)

   dfs = [pd.read_parquet(f) for f in FILI]
   trange = [df.index.max() - df.index[0] for df in dfs]
   print(trange)
   totyrs = [tr.total_seconds() / 86400 / 365 for tr in trange]
   print(totyrs)
   #exit()

   bininfo = make_bin_grid_1d(latbnds=(-90., 90.), dlat=DLAT)

   plt.rc('font', size=14)

   uniq = plot_lat_binned(dfs, bininfo, totyrs, unique_tracks_of_storm_1d, 'unique storms', 'unique_track_dens_%.1f.png' % DLAT)
   h6all = plot_lat_binned(dfs, bininfo, totyrs, sixhrly_fixes_of_storm_1d, '6-hourly fixes', '6hr_track_dens_%.1f.png' % DLAT)
   h6hurr = plot_lat_binned(dfs, bininfo, totyrs, sixhrly_fixes_of_storm_1d, '6-hourly hurricane fixes', '6hr_track_dens_hurr_%.1f.png' % DLAT, wspd_thresh=33)
   h6maj = plot_lat_binned(dfs, bininfo, totyrs, sixhrly_fixes_of_storm_1d, '6-hourly major hurricane fixes', '6hr_track_dens_major_%.1f.png' % DLAT, wspd_thresh=50)
   ace = plot_lat_binned(dfs, bininfo, totyrs, bin_ace_of_storm_1d, 'ACE [$10^4$ kt$^2$] ', 'ace_binned_%.1f.png' % DLAT)

   outds = xr.Dataset(data_vars=dict(uniq=uniq, h6all=h6all, h6hurr=h6hurr, h6maj=h6maj, ace=ace), attrs=dict(dlat=DLAT))
   outds.to_netcdf(os.path.join(DOUT, 'tcdens.nc'))

def plot_lat_binned(dfs, bininfo, totyrs, varfunc, ylabel, filo, norm='per million km2 per yr', **fkwargs):
   pltdat = [np.zeros_like(bininfo[1]) for _ in dfs]
   for ii, df in enumerate(dfs):
      for stm in df['stmnum'].unique():
         pltdat[ii] += varfunc(df[df['stmnum'] == stm], bininfo, **fkwargs)
      print(ylabel, '\t', pltdat[ii].sum())
      pltdat[ii] = pltdat[ii] / (bininfo[2] / 1e6 / 1e3**2) / totyrs[ii]
      
      plt.plot(np.sin(np.deg2rad(bininfo[1])), pltdat[ii], label=labels[ii])
      plt.xticks(np.sin(np.deg2rad(LATLABS)), labels=LATLABS.astype(np.int_))
      plt.xlabel('lat')
      if varfunc == sixhrly_fixes_of_storm_1d:
         plt.ylim(*FIXYLIM)
      plt.ylabel('%s %s' % (ylabel, norm))
   plt.legend(fontsize=12)
   plt.savefig(os.path.join(DOUT, filo), bbox_inches='tight')
   plt.close()

   return xr.DataArray(np.array(pltdat), dims=['run', 'lat'], coords=dict(run=labels, lat=bininfo[1]), name=ylabel, attrs=dict(norm=norm))

#make an array of lat bin interfaces
#In: latitude bounds [deg], desired lat spacing [deg] at the equator
def make_bin_grid_1d(latbnds=(-50., 50.), dlat=0.5):
   sinbnds = (np.sin(np.deg2rad(latbnds[0])), np.sin(np.deg2rad(latbnds[1])))
   dsin = np.sin(np.deg2rad(dlat))
   true_dsin = (sinbnds[1] - sinbnds[0]) / ((sinbnds[1] - sinbnds[0]) // dsin + 1) #choose a spacing so that the interval is evenly divided

   bini = np.arange(sinbnds[0], sinbnds[1] + 1e-3, true_dsin) #bin interfaces in sin(lat) space
   binm = (bini[1:] + bini[:-1]) / 2
   area = 2 * np.pi * a_e**2 * true_dsin
   return np.rad2deg(np.arcsin(bini)), np.rad2deg(np.arcsin(binm)), area

#Binary search latitude bins and accumulate the desired quantity
#In: dataframe subset for 1 storm, bininfo (bin interfaces, bin midpoints, bin area [m^2]), column name, func to apply to column values, accum data type, return data type
def accum_bin_map_1d(stmdf, bininfo, col, afunc, dtype=np.float64, rettype=np.float64):
   outarr = np.zeros_like(bininfo[1], dtype=dtype)
   accvals = stmdf[col].apply(afunc)
   for ii in stmdf.index:
      if np.isnan(stmdf['lat'][ii]):
         continue
      iidx = np.searchsorted(bininfo[0], stmdf['lat'][ii])
      if iidx == 0 or iidx == bininfo[0].shape[0] - 1:
         continue
      outarr[iidx - 1] += accvals[ii]
   return outarr.astype(rettype)

#bool-to-int array of whether the storm hit the bin with midpoint
#i.e., return should be 0 or 1 divided by bin area
#In: dataframe subset for 1 storm, bininfo (bin interfaces, bin midpoints, bin area [m^2])
def unique_tracks_of_storm_1d(stmdf, bininfo):
   return accum_bin_map_1d(stmdf, bininfo, 'lat', lambda x: 1, dtype=np.bool_, rettype=np.int_)

def sixhrly_fixes_of_storm_1d(stmdf, bininfo, wspd_thresh=0):
   return accum_bin_map_1d(stmdf, bininfo, 'wspd', lambda wspd: int(wspd >= wspd_thresh), dtype=np.int_, rettype=np.int_)

def bin_ace_of_storm_1d(stmdf, bininfo):
   ace = lambda wspd: 1e-4 * (MS2KT * wspd)**2
   return accum_bin_map_1d(stmdf, bininfo, 'wspd', ace)

if __name__ == '__main__':
   main()
