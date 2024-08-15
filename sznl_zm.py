#Joshua Pan orig from Oct 2023 for docn runs
#plot and compare zonal-mean time-averaged quantities

#TOPDIR = '/glade/scratch/zarzycki/archive/'
#DIR1 = 'QPC5-ne30np4-aquap10-seed3x3/atm/hist'
#DIR2 = 'QPC5-ne30np4-aquap10-unseed/atm/hist'
#FN = 'QPC5-ne30np4-aquap10-*.cam.h0.*regrid.nc'
OUTDIR = 'linevslat_0812/'
HISTDIMS = set(['time', 'lat', 'lon']) #cam
#HISTDIMS = set(['time', 'xh', 'yh']) #mom6 hm

import proj3
import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pltsettings

#print(dir(proj3))
ARCHV = proj3.ARCHV
CASE = proj3.CASE
H0 = proj3.H0
HISTS = proj3.HISTS
tdel = proj3.tdel
wtm = proj3.weighted_temporal_mean

LATLAB = np.array([-90., -60., -30., 0., 30., 60., 90.])
lncolors = plt.cm.jet(np.linspace(0, 1, 12))

def main():
   pltsettings.set1()

   print('Opening datasets...')
   ds1 = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H0))#!OCEAN .rename(dict(xh='lon', yh='lat'))
   ds1 = ds1.assign_coords(coords=dict(time=ds1.time - tdel(days=1))) #timestamp of e.g., 0001-01.nc is Feb 1, so subtract a month from time coord
   #ds1 = ds1.assign(dict(U200=ds1.U.sel(lev=200, method='nearest')))
   #!OCEAN ds1 = ds1.isel(z_l=0)
   print(ds1)

   for dv in ds1.data_vars:
      #if str(dv) != 'TREFHT':
      #   print('Skipping %s...' % dv)
      #   continue
      if set(ds1[dv].dims) == HISTDIMS:
         print('Plotting %s...' % dv)
         monmeans = ds1[dv].groupby('time.month').mean()
         lines = monmeans.mean(dim='lon')
         sinlat = np.sin(np.deg2rad(lines.lat))
         for mo in lines.month:
            plt.plot(sinlat, lines.sel(month=mo), label=mo.values, color=lncolors[mo.values-1])
         #plt.plot(sinlat, line1.values)
         if str(dv) == 'TS':
            plt.hlines(273.15 + 26.5, -1, 1, colors='black', linestyles='dashed')
         plt.xlabel('Lat [°]')
         plt.ylabel(dv)
         plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=8))
         plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB)
         plt.savefig(os.path.join(OUTDIR, '%s_sznl.png' % dv), bbox_inches='tight')
         #plt.tight_layout()
         #plt.show()
         plt.close()
      else:
         print('Skipping %s...' % dv)

   print('%s done.' % sys.argv[0])

if __name__ == '__main__':
   main()
