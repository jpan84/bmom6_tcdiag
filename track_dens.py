#Joshua Pan Mar 2025
#Uses parquet output of traj_stats_bmom_SE.py
#CL arg: parquet file

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

FILI = sys.argv[1]

SSHS = np.array([0, 17, 33, 43, 49, 58, 70, 9999])
COLS = [None, 'ncol', 'lon', 'lat', 'pres', 'wspd', 'year', 'month', 'day', 'hour']
TYPS = [None, int, float, float, float, float, int, int, int, int]
CTYP = {COLS[i+1]: TYPS[i+1] for i in range(len(COLS[1:]))}

XLIMS = dict(pmins=(8.5e4,1.01e5), ace=(0,80), maxu=(15,120), genlon=(0,360), genlat=(-40, 40))
YLIMS = dict(pmins=(0, 6e-2))
clabelkwargs = {'inline': 1, 'fontsize': 10, 'colors': 'black', 'fmt': '%.1e'}

def main():
   extidx = FILI.rfind('.')
   DOUT = FILI[:extidx] + '_density'
   if not os.path.exists(DOUT):
      os.makedirs(DOUT)

   df = pd.read_parquet(FILI)

   bininfo = make_bin_grid_1d()

   for stm in df['stmnum'].unique():
      tr_dens = track_density_of_storm_1d(df[df['stmnum'] == stm], bininfo)
      ace_stm = bin_ace_of_storm_1d(df[df['stmnum'] == stm], bininfo)

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

#track density in # of storms per 1e6 km**2
#i.e., return should be 0 or 1 divided by bin area
#In: dataframe subset for 1 storm, bininfo (bin interfaces, bin midpoints, bin area [m^2])
def track_density_of_storm_1d():


if __name__ == '__main__':
   main()
