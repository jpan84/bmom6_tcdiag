from paths import CASENAMES, IXHORS, CTLIX
import pandas as pd
import numpy as np
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt

TRAJPQS = ['TC_preprocess/trajectories.txt.%s.parquet' % lc for lc in CASENAMES] [1:]
# Grid resolution for KDE evaluation

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

   plt.rcParams['figure.figsize'] = (24, 10)
   fig, axes = plt.subplots(2, 3, sharex=True, sharey=True)

   for ii, ax in enumerate(axes.ravel()):
      if ii in [3, 4]: continue
      ixh = IXHORS[ii] - 1
      ctlix = CTLIX - 1
      isctl = ixh == ctlix
      pltar = binout[ixh][0]
      if not isctl:
         pltar -= binout[ctlix][0]
      #ax.pcolormesh(DOY_GC, LAT_GC, binout[ixh][0].T, shading='nearest')
      csf = ax.contourf(DOY_GC, LAT_GC, pltar.T, cmap='inferno' if isctl else 'bwr')
      cb = plt.colorbar(csf, ax=ax)

   plt.show()

if __name__ == '__main__':
   main()
