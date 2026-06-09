import os
from paths import CASENAMES, IXHORS, CTLIX, ARCHRT
import cftime
import pandas as pd
import numpy as np
import xarray as xr
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt

TRAJPQS = ['TC_preprocess/trajectories.txt.%s.parquet' % lc for lc in CASENAMES] [1:]
SSTZM = 'atm/uxzm_h1i_noncons_-20.0_20.0_0.25_SST.nc'

DOY_GE = np.arange(0.5, 366, 1)
MU_GE = np.linspace(-1, 1, 200) #mu=sinlat

DOY_GE = np.arange(0.5, 366, 7)
MU_GE = np.linspace(-1, 1, 40) #mu=sinlat

MU2LAT = lambda x: np.rad2deg(np.arcsin(np.clip(x, -1., 1.)))
LAT2MU = lambda x: np.sin(np.deg2rad(np.clip(x, -90., 90.)))
GE2GC = lambda x: (x[1:] + x[:-1]) / 2 #grid edges to midpoints

LAT_GE = MU2LAT(MU_GE)
LAT_GC = GE2GC(LAT_GE)
DOY_GC = GE2GC(DOY_GE)
MU_GC = GE2GC(MU_GE)

TICKDATES = [cftime.datetime(1, ii, 1, calendar='noleap') for ii in range(1, 13)]
TICKDOY = [int(dt.strftime('%-j')) for dt in TICKDATES]

def main():
   dfs = [pd.read_parquet(tp) for tp in TRAJPQS]
   for df in dfs:
      if type(df.index) != pd.core.indexes.datetimes.DatetimeIndex:
         df['dt'] = pd.to_datetime(df['dt'])
         df.set_index('dt', inplace=True)
   trange = [df.index.max() - df.index[0] for df in dfs]
   totyrs = [tr.total_seconds() / 86400 / 365 for tr in trange]
   print(totyrs)

   gen = [df[df['isgen']] for df in dfs]
   print(gen[0])

   binout = [binned_statistic_2d(gdf['doy'], LAT2MU(gdf['lat']), gdf['lat'],\
             statistic='count', bins=[DOY_GE, MU_GE])\
             for gdf in gen]
   #genct, _, _, binix

   sstdss = []
   for ar in ARCHRT:
      #print(os.path.join(ar, SSTZM))
      try:
         sstdss.append(xr.open_dataset(os.path.join(ar, SSTZM)))
      except:
         sstdss.append(None)
   #print(sstdss)

   maxes = [None if ds is None else ds['SST'].idxmax('latitudes') for ds in sstdss]
   del maxes[0]
   #print(maxes)
   #exit()
   #sgn_chg = np.sign(maxes).diff(dim='time')
   #chg_idx = sgn_chg['time'].where(sgn_chg != 0, drop=True)
   #half_year = (chg_idx.dt.month > 6).astype(int)
   #grp_key = chg_idx.dt.year + half_year * 0.5
   #chg_idx.coords["half_years"] = grp_key
   ##print(chg_idx)
   #first_flip = chg_idx.groupby('half_years').first()
   ##print(first_flip)
   #flipdoys = first_flip.dt.dayofyear
   #sprg, wntr = flipdoys.where(flipdoys > 182, drop=True).mean(), flipdoys.where(flipdoys <= 182, drop=True).mean()

   plt.rcParams['figure.figsize'] = (24, 10)
   fig, axes = plt.subplots(2, 3, sharex=True, sharey=True)

   for ii, ax in enumerate(axes.ravel()):
      if ii in [3, 4]: continue
      ixh = IXHORS[ii] - 1
      ctlix = CTLIX - 1
      isctl = ixh == ctlix
      pltar = binout[ixh][0]

      pltts = maxes[ixh] #None if maxes[ixh] is None else maxes[ixh]['SST'].idxmax('latitudes')
      sp, wn = av_SST_flip_doy(pltts)

      clevs = np.arange(0, 9) if isctl else np.arange(-8, 9)
      if not isctl:
         pltar -= binout[ctlix][0]
      #ax.pcolormesh(DOY_GC, LAT_GC, binout[ixh][0].T, shading='nearest')
      csf = ax.contourf(DOY_GC, LAT_GC, pltar.T, cmap='inferno' if isctl else 'bwr', levels=clevs)
      cb = plt.colorbar(csf, ax=ax)
      ax.scatter(pltts['time'].dt.dayofyear, pltts, alpha=0.2, c='lime', s=0.5) if pltts is not None else None
      [ax.axvline(fd, c='gray', linestyle='dashed') if fd is not None else None for fd in [sp, wn]]

      ax.set_ylim(-30, 30)
      ax.set_xticks(TICKDOY, [dt.strftime('%m-%d') for dt in TICKDATES], rotation=45)
      ax.tick_params(axis='both', labelleft=True, right=True, top=True)

   plt.show()

def av_SST_flip_doy(maxlat):
   if maxlat is None: return None, None

   sgn_chg = np.sign(maxlat).diff(dim='time')
   chg_idx = sgn_chg['time'].where(sgn_chg != 0, drop=True)
   half_year = (chg_idx.dt.month > 6).astype(int)
   grp_key = chg_idx.dt.year + half_year * 0.5 
   chg_idx.coords["half_years"] = grp_key
   #print(chg_idx)
   first_flip = chg_idx.groupby('half_years').first()
   #print(first_flip)
   flipdoys = first_flip.dt.dayofyear
   sprg, wntr = flipdoys.where(flipdoys > 182, drop=True).mean(), flipdoys.where(flipdoys <= 182, drop=True).mean()
   return sprg, wntr

if __name__ == '__main__':
   main()
