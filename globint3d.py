#Joshua Pan Oct 2023
#compute globally area-weighted averaged quantities from unstructured grid

ARCHV = '/glade/derecho/scratch/jpan/archive/'
HISTS = 'atm/hist/'
H0 = r'*h0a.[0-9]*.nc'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl'
GRIDFN = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
OUTDIR = './globavtraces_250417_ctrl'

CLIPMO = 0
tunits = 'common_years since 0000-01-01'

import os
import glob
import cftime
import numpy as np
import uxarray as ux
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#import pltsettings
import consts as c

a = c.a
GA = 4 * np.pi * a**2

RAW = ['Q', 'DCQ', 'T', 'U', 'V']
#UDEF = ['KE', 'GP', 'SH', 'LH', 'DIAB',
UDEF = dict(KE=[(1, 'UU'), (1, 'VV'), 'sum'], )

plt.rc('font', size=16)

def main():
   if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

   pt = os.path.join(ARCHV, CASE, HISTS, H0)
   ds = ux.open_mfdataset(os.path.join(GRIDDIR, GRIDFN), pt)
   flt = cftime.date2num(ds.time, tunits)
   print(flt)


   gav, units = None, None
   for var in vrs:
      print('Working on', var)
      if var == 'PRECT':
         gav = globav(ds, 'PRECC') + globav(ds, 'PRECL')
         units = ds['PRECC'].units
      elif var == 'RESTOM':
         gav = globav(ds, 'FSNT') - globav(ds, 'FLNT')
         units = ds['FSNT'].units
      elif var == 'NCF':
         gav = globav(ds, 'LWCF') + globav(ds, 'SWCF')
         units = ds['LWCF'].units
      else:
         gav = globav(ds, var)
         units = ds[var].units
      with open(os.path.join(OUTDIR, 'globavgs.txt'), 'a') as f:
         if var == vrs[0]:
            print('\n', pt, file=f)
            print(gav.isel(time=slice(CLIPMO, None)).shape, '\n', file=f)
         #print(var, tav(gav.isel(time=slice(CLIPMO,))), file=f)
         print(var, weighted_temporal_mean(gav.isel(time=slice(CLIPMO, None))), file=f)
         f.close()
      #print(gav.time)
      #print(weighted_temporal_mean(gav))
      plt.plot(flt, gav.values)
      #plt.legend()
      plt.title(var)
      plt.xlabel(tunits)# [%s]' % str(flt.units))
      plt.ylabel(units)
      plt.savefig(os.path.join(OUTDIR, '%s.png' % var), bbox_inches='tight')
      plt.close()

   print('globavg.py done.')

#from mf_dataset with uxgrid, compute global average of var
def globav(ds, varname):
   v = ds[varname]
   #if 'lev' in v.dims:
   #   print(v.lev.isel(lev=-1))
   #   v = v.isel(lev=-1)
   vxarea = v * ds.uxgrid.face_areas * a**2 #m2 * [v], area times var at each face and time
   gav = vxarea.sum(dim='n_face') / (4 * np.pi * a**2)
   return gav

#compute time avg of global avg
#assume uniform timestep
def tav(gav):
   dt = np.repeat(gav.time.values[1] - gav.time.values[0], gav.time.shape[0])
   print(gav.shape, dt.shape)
   tchunks = gav * dt #s * [v], time integral of global avg var each time interval
   return (tchunks.sum() / dt.sum()).values

#adapted from https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
#weight monthly avgs to correctly compute annual avgs
def weighted_temporal_mean(da):
   """
   weight by days in each month
   """
   # Determine the month length
   month_length = da.time.dt.days_in_month
   #print(month_length)

   '''
   # Calculate the weights
   wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
   print(wgts)

   # Make sure the weights in each year add up to 1
   np.testing.assert_allclose(wgts.groupby("time.year").sum(xr.ALL_DIMS), 1.0)

   # Subset our dataset for our variable
   #obs = ds[var]
   obs = da

   # Setup our masking for nan values
   cond = obs.isnull()
   ones = xr.where(cond, 0.0, 1.0)

   # Calculate the numerator
   obs_sum = (obs * wgts).resample(time="AS").sum(dim="time")

   # Calculate the denominator
   ones_out = (ones * wgts).resample(time="AS").sum(dim="time")

   # Return the weighted average
   return obs_sum / ones_out
   '''

   num = (da * month_length).sum(dim='time')
   den = month_length.sum(dim='time')
   return num / den

if __name__ == '__main__':
   main()
