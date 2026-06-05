from paths import CASENAMES, IXHORS
import pandas as pd
import numpy as np
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt

TRAJPQS = ['TC_preprocess/trajectories.txt.%s.parquet' % lc for lc in CASENAMES] [1:]
# Grid resolution for KDE evaluation

DOY_GE = np.arange(0.5, 366, 1)
MU_GE = np.linspace(-1, 1, 200) #mu=sinlat

MU2LAT = lambda x: np.rad2deg(np.arcsin(np.clip(x, -1., 1.)))
LAT2MU = lambda x: np.sin(np.deg2rad(np.clip(x, -90., 90.)))
GE2GC = lambda x: (x[1:] + x[:-1]) / 2 #grid edges to midpoints

LAT_GE = MU2LAT(MU_GE)
LAT_GC = GE2GC(LAT_GE)
DOY_GC = GE2GC(DOY_GE)
MU_GC = GE2GC(MU_GE)


def main():
   dfs = [pd.read_parquet(tp) for tp in TRAJPQS]
   gen = [df[df['isgen']] for df in dfs]
   print(gen[0])

   binout = [binned_statistic_2d(gdf['doy'], LAT2MU(gdf['lat']), gdf['lat'],\
             statistic='count', bins=[DOY_GE, MU_GE])\
             for gdf in gen]
   #genct, _, _, binix

   fig, axes = plt.subplots(2, 3, sharex=True, sharey=True)

   for ii, ax in enumerate(axes.ravel()):
      if ii in [3, 4]: continue
      ixh = IXHORS[ii] - 1
      print(binout[ixh][0].max())
      ax.pcolormesh(DOY_GC, LAT_GC, binout[ixh][0].T, shading='nearest')

   plt.show()

if __name__ == '__main__':
   main()
