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

#FILI = sys.argv[1]
FILI = ['8e-5.parquet', '250206_1degeqm885.parquet', '250321_seed1x1.parquet']
labels = ['unseed', 'ctrl', 'seed']
DOUT = '250402_density'

def main():
   #extidx = FILI.rfind('.')
   #DOUT = FILI[:extidx] + '_density'
   if not os.path.exists(DOUT):
      os.makedirs(DOUT)

   dfs = [pd.read_parquet(f) for f in FILI]
   trange = [df.index.max() - df.index[0] for df in dfs]
   totyrs = [tr.total_seconds() / 86400 / 365 for tr in trange]

   bininfo = make_bin_grid_1d()
   tr_dens_1d = [np.zeros_like(bininfo[1]) for _ in dfs]
   bin_ace_1d = [np.zeros_like(bininfo[1]) for _ in dfs]

   plt.rc('font', size=16)
   for ii, df in enumerate(dfs):
      for stm in df['stmnum'].unique():
         #print('Storm #%d' % stm)
         tr_dens_1d[ii] += track_density_of_storm_1d(df[df['stmnum'] == stm], bininfo)   

      tr_dens_1d[ii] /= (bininfo[2] / 1e6 / 1e3**2) / totyrs[ii] #track density in # of storms per 1e6 km**2 per year

      plt.plot(bininfo[1], tr_dens_1d[ii], label=labels[ii])
      plt.xlabel('lat')
      plt.ylabel('Storms per million km2 per yr')
   plt.legend(fontsize=12)
   plt.savefig(os.path.join(DOUT, 'track_dens.png'), bbox_inches='tight')
   plt.close()

   for ii, df in enumerate(dfs):
      for stm in df['stmnum'].unique():
         bin_ace_1d += bin_ace_of_storm_1d(df[df['stmnum'] == stm], bininfo)

      bin_ace_1d[ii] /= (bininfo[2] / 1e6 / 1e3**2) / totyrs[ii]

      plt.plot(bininfo[1], bin_ace_1d[ii], label=labels[ii])
      plt.xlabel('lat')
      plt.ylabel('ACE [$10^{-4}$ kt$^2$] per million km2 per yr')
   plt.legend(fontsize=12)
   plt.savefig(os.path.join(DOUT, 'ace_binned.png'), bbox_inches='tight')
   plt.close()

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

#bool-to-int array of whether the storm hit the bin with midpoint
#i.e., return should be 0 or 1 divided by bin area
#In: dataframe subset for 1 storm, bininfo (bin interfaces, bin midpoints, bin area [m^2])
def track_density_of_storm_1d(stmdf, bininfo):
   boolhit = np.full_like(bininfo[1], False)
   for lat in stmdf['lat']:
      #print(lat)
      iidx = np.searchsorted(bininfo[0], lat) #edge case of lat==upper bound?
      boolhit[iidx - 1] = True
   return boolhit.astype(np.int_)

def bin_ace_of_storm_1d(stmdf, bininfo):
   binned_ace = np.zeros_like(bininfo[1])
   ace = 1e-4 * (MS2KT * stmdf['wspd'])**2
   for ii in stmdf.index:
      #print(lat)
      iidx = np.searchsorted(bininfo[0], stmdf['lat'][ii]) #edge case of lat==upper bound?
      binned_ace[iidx - 1] += ace[ii]
   return binned_ace

if __name__ == '__main__':
   main()
