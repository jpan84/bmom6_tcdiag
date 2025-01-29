#Joshua Pan Oct 2023
#compute globally area-weighted averaged quantities from unstructured grid

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.ne120pg3_sx0.66av1.aqua.zmtau.0826/'
HISTS = 'atm/hist/'
H0 = r'*h0a.[0-9]*.nc'
GRIDDIR = '/glade/p/cesmdata/inputdata/share/scripgrids'
GRIDFN = 'ne120pg3_scrip_170628.nc'
OUTDIR = './globavtraces_0826_zmtau'


CASE = 'b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.1121_yr1098/'
HISTS = 'atm/hist/'
H0 = r'*h0a.[0-9]*.nc'
GRIDFN = 'ne30np4_091226_pentagons.nc'
OUTDIR = './globavtraces_1121_yr1098'

CASE = 'b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.250128_h80l894/'
HISTS = 'atm/hist/'
H0 = r'*h0a.[0-9]*.nc'
GRIDFN = 'ne30np4_091226_pentagons.nc'
OUTDIR = './globavtraces_250128_h80l894'

import os
import glob
import numpy as np
import uxarray as ux
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#import pltsettings
import consts as c

a = c.a

vrs = ['TS', 'FSNT', 'FLNT', 'RESTOM', 'PRECC', 'PRECL', 'PRECT', 'QFLX', 'PS', 'TMQ', 'LWCF', 'SWCF', 'CLDTOT']
plt.rc('font', size=16)

def main():
   pt = os.path.join(ARCHV, CASE, HISTS, H0)
   #print(os.path.join(GRIDDIR, GRIDFN))
   #print(pt)
   #exit()
   ds = ux.open_mfdataset(os.path.join(GRIDDIR, GRIDFN), pt)
   print(ds.time) 

   '''
   #pltsettings.set1()
   args = [(TOPDIR, DIR1, FN), (TOPDIR, DIR2, FN), (TOPDIR3, DIR3, FN)]
   aliases = ['seed', 'unseed', 'ctrl']
   paths = [sorted(glob.glob(os.path.join(*a))) for a in args]
   dss = [ux.open_mfdataset(os.path.join(GRIDDIR, GRIDFN), ps) for ps in paths]
   '''

   gav, units = None, None
   for var in vrs:
      if var == 'PRECT':
         gav = globav(ds, 'PRECC') + globav(ds, 'PRECL')
         units = ds['PRECC'].units
      elif var == 'RESTOM':
         gav = globav(ds, 'FSNT') - globav(ds, 'FLNT')
         units = ds['FSNT'].units
      else:
         gav = globav(ds, var)
         units = ds[var].units
      with open(os.path.join(OUTDIR, 'globavgs.txt'), 'a') as f:
         if var == vrs[0]:
            print(pt, '\n', file=f)
         print(var, tav(gav), file=f)
         f.close()
      #print(gav.time)
      #print(weighted_temporal_mean(gav))
      plt.plot(gav.time, gav.values)
      plt.legend()
      plt.title(var)
      plt.xlabel('Time [%s]' % str(gav.time.units))
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

   # Calculate the weights
   wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()

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

if __name__ == '__main__':
   main()
