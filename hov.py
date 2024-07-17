#Joshua Pan

import proj3
import os
import sys
import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pltsettings

#print(dir(proj3))
ARCHV = proj3.ARCHV
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0711dlyout'
H1 = r'*h1.[0-9]*.nc'
HISTS = proj3.HISTS
tdel = proj3.tdel

LATBNDS = (-10, 10)
LONBNDS = (0, 90)
TBNDS = None

def main():
   pltsettings.set1()

   print('Opening datasets...')
   ds1 = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H1))
   TBNDS = (ds1.time[0], ds1.time[0] + tdel(days=365))

   pltvar = ds1.PRECT
   with ProgressBar():
      pltvar = pltvar.sel(lat=slice(*LATBNDS), lon=slice(*LONBNDS), time=slice(*TBNDS)).mean(dim='lat').load()

   [print(tt) for tt in pltvar.time]
   [print(tt - pltvar.time[0]).values for tt in pltvar.time.values]
   #[print(type(tt - pltvar.time[0])) for tt in pltvar.time.values]
   taxis = [(tt - pltvar.time[0]).values.days for tt in pltvar.time.values]

   plt.contourf(pltvar.lon, [(tt - pltvar.time[0]).values.days for tt in pltvar.time.values], pltvar.values)
   plt.colorbar()
   plt.savefig('hov_test.png')
   plt.show()
   plt.close()


   exit()
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
