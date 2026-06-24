import os
from paths import CASENAMES, IXHORS, CTLIX, ARCHRT, ALIA
import cftime
import pandas as pd
import numpy as np
import xarray as xr
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt
import matplotlib.colors as colors

TRAJPQS = ['TC_preprocess/trajectories.txt.%s.parquet' % lc for lc in CASENAMES]# [1:]
SSTZM = 'atm/uxzm_h1i_noncons_-18.0_18.0_0.1_SST.nc'

DOY_GE = np.arange(0.5, 366, 1)
MU_GE = np.linspace(-1, 1, 200) #mu=sinlat

DOY_GE = np.arange(0.5, 366, 5)
MU_GE = np.linspace(-1, 1, 80) #mu=sinlat

MU2LAT = lambda x: np.rad2deg(np.arcsin(np.clip(x, -1., 1.)))
LAT2MU = lambda x: np.sin(np.deg2rad(np.clip(x, -90., 90.)))
GE2GC = lambda x: (x[1:] + x[:-1]) / 2 #grid edges to midpoints
#GE2GD = lambda x: (x[1:] - x[:-1]) #grid edges to grid deltas

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

   bandarea = 2 * np.pi * 6.371e6**2 / 1e6 / 1e3**2 * np.diff(MU_GE)
   dday = np.diff(DOY_GE)

   binout = [binned_statistic_2d(gdf['doy'], LAT2MU(gdf['lat']), gdf['lat'],\
             statistic='count', bins=[DOY_GE, MU_GE])\
             for gdf in gen]
   binnormed = [bo[0] / totyrs[ii] / dday[:, None] / bandarea[None, :] for ii, bo in enumerate(binout)]
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

   plt.rcParams['figure.figsize'] = (24, 10)
   plt.rc('font', size=16)
   fig, axes = plt.subplots(2, 3, sharex=True, sharey=True)

   for ii, ax in enumerate(axes.ravel()):
      if ii in [4]:
         ax.set_axis_off()
         continue
      ixh = IXHORS[ii] #- 1
      ctlix = CTLIX #- 1
      isctl = ixh == ctlix
      pltar = binnormed[ixh]

      pltts = maxes[ixh] #None if maxes[ixh] is None else maxes[ixh]['SST'].idxmax('latitudes')
      sp, wn = av_SST_flip_doy(pltts)

      clevs = np.arange(0, .031, .003) if isctl else np.arange(-.03, .031, .003)
      #bnorm = colors.BoundaryNorm(clevs, 256)
      
      if not isctl:
         pltar -= binnormed[ctlix]#[0]
      #ax.pcolormesh(DOY_GC, LAT_GC, binout[ixh][0].T, shading='nearest')

      # 1. Grab your base colormap resampled to your exact number of levels
      nlevs = len(clevs) - 1
      cmap = plt.get_cmap('inferno' if isctl else 'seismic').resampled(nlevs)

      if not isctl:
         # 2. Extract the RGBA array from the colormap
         rgba_colors = cmap(np.arange(nlevs))
         
         # 3. Dynamically find which bins contain zero based on your 'clevs'
         # This finds any bin where the left edge is <= 0 AND the right edge is >= 0
         zero_bins = np.where((clevs[:-1] <= 0) & (clevs[1:] >= 0))[0]
         # Change line 98 to check if either boundary touches zero:
         zero_bins = np.where((clevs[:-1] == 0) | (clevs[1:] == 0))[0]
         # Use np.isclose to catch values that are effectively 0, bypassing float precision errors
         zero_bins = np.where(np.isclose(clevs[:-1], 0, atol=1e-7) | np.isclose(clevs[1:], 0, atol=1e-7))[0]
         
         # 4. Turn exactly those center bins pure white
         rgba_colors[zero_bins] = [1, 1, 1, 1]
         cmap = colors.ListedColormap(rgba_colors)
      bnorm = colors.BoundaryNorm(clevs, ncolors=nlevs)

      csf = ax.pcolormesh(DOY_GE, MU_GE, pltar.T, shading='flat', cmap=cmap, norm=bnorm)#, levels=clevs)
      cb = plt.colorbar(csf, ax=ax)
      ax.scatter(pltts['time'].dt.dayofyear, LAT2MU(pltts), alpha=0.2, c='lime', s=0.5) if pltts is not None else None
      [ax.axvline(fd, c='C1', linestyle='dashed') if fd is not None else None for fd in [sp, wn]]

      ax.set_xlim(0.5, 365.5)
      ax.set_ylim(-0.5, 0.5)
      ax.set_xticks(TICKDOY, [dt.strftime('%m-%d') for dt in TICKDATES], rotation=45)
      ax.set_yticks(LAT2MU(np.arange(-30, 31, 10)), np.arange(-30, 31, 10))
      ax.tick_params(axis='both', labelleft=True, labelbottom=True, right=True, top=True)
      ax.set_title(ALIA[ixh])
      ax.set_title('(%s)' % chr(ord('a') + ii))

      print(sp, wn)
      fddt = [cftime.num2date(fd, units='days since 0000-12-31', calendar='noleap') if not np.isnan(fd) else None for fd in [sp, wn]]
      fdmmdd = [ft.strftime('%m-%d') if ft is not None else None for ft in fddt]
      [ax.text(fd - 50, -.1, fdmmdd[jj], c='C1') if not np.isnan(fd) else None for jj, fd in enumerate([sp, wn])]

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
