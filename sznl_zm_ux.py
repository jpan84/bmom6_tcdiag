#Joshua Pan orig from Oct 2023 for docn runs
#Updated May 2025 to use uxarray non-conservative zonal mean
#plot and compare zonal-mean time-averaged quantities

#TOPDIR = '/glade/scratch/zarzycki/archive/'
#DIR1 = 'QPC5-ne30np4-aquap10-seed3x3/atm/hist'
#DIR2 = 'QPC5-ne30np4-aquap10-unseed/atm/hist'
#FN = 'QPC5-ne30np4-aquap10-*.cam.h0.*regrid.nc'
OUTDIR = 'linevslat_250416_seed1x1/'
HISTDIMS = set(['time', 'n_face']) #cam-SE

import os
import sys
import cftime
import datetime as dt
import numpy as np
import uxarray as ux
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pltsettings

### hist file params
ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1'
HISTS = 'atm/hist/'
H0 = '*.cam.h0a.*.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

grpby = 'month' #'season'
zmlats = (-90, 90, 0.5)
LATLAB = np.array([-90., -60., -30., 0., 30., 60., 90.])
lncolors = plt.cm.jet(np.linspace(0, 1, 12 if grpby == 'month' else 4 if grpby == 'season'))
#TODO: allow diffing between cases and selecting of months/seasons

def main():
   if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

   pltsettings.set1()

   print('Opening datasets...')
   
   ds1 = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HISTS, H0))
   #ds1 = ds1.assign_coords(coords=dict(time=ds1.time - dt.timedelta(days=1))) #timestamp of e.g., 0001-01.nc is Feb 1, so subtract a month from time coord
   #ds1 = ds1.assign(dict(U200=ds1.U.sel(lev=200, method='nearest')))
   #!OCEAN ds1 = ds1.isel(z_l=0)
   print(ds1)

   for dv in ds1.data_vars:
      #if str(dv) != 'TREFHT':
      #   print('Skipping %s...' % dv)
      #   continue
      if set(ds1[dv].dims) == HISTDIMS:
         print('Plotting %s...' % dv)
         tmeans = None
         if grpby == 'month':
            tmeans = ds1[dv].groupby('time.%s' % grpby).mean()
            tmeans = tmeans.assign_coords(month=np.arange(12) + 1)
         elif grpby == 'season':
            tmeans = season_mean(ds1[dv])
         lines = tmeans.zonal_mean(lat=zmlats)
         sinlat = np.sin(np.deg2rad(lines.latitudes))
         print(lines)
         print(lines.coords)
         print(lines.month)
         for tt in lines.coords[grpby]:
            plt.plot(sinlat, lines.sel(**{grpby: tt}), label=tt.values, color=lncolors[tt.values])
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

#https://docs.xarray.dev/en/stable/examples/monthly-means.html
#Seasonal means from monthly means, properly weighted by month length
def season_mean(ds, calendar="noleap"):
   # Make a DataArray with the number of days in each month, size = len(time)
   month_length = ds.time.dt.days_in_month

   # Calculate the weights by grouping by 'time.season'
   weights = (
       month_length.groupby("time.season") / month_length.groupby("time.season").sum()
   )

   # Test that the sum of the weights for each season is 1.0
   np.testing.assert_allclose(weights.groupby("time.season").sum().values, np.ones(4))

   # Calculate the weighted average
   return (ds * weights).groupby("time.season").sum(dim="time")

if __name__ == '__main__':
   main()
