#Joshua Pan Oct 2023
#compute globally area-weighted averaged quantities from unstructured grid

ARCHV = '/glade/derecho/scratch/jpan/archive/'
HISTS = 'atm/hist/'
H0 = r'*h0a.[0-9]*.nc'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl'
GRIDFN = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
OUTDIR = './globinteg_250417_ctrl'

CLIPMO = 0
tunits = 'common_years since 0000-01-01'

import os
import glob
import cftime
import numpy as np
import uxarray as ux
#import xarray as xr
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#import pltsettings
import consts as c

a = 6.371e6
GA = 4 * np.pi * a**2
P0 = 1e5
g = 9.81

RAW = ['Q', 'DCQ', 'T', 'U', 'V']
#UDEF = ['KE', 'GP', 'SH', 'LH', 'DIAB',
UDEF = dict(KE=[(1, 'UU'), (1, 'VV'), 'sum'], GP=[(g, 'Z3'), 'sum'], MASS=[(1, ), 'sum'])
VARS = RAW + list(UDEF.keys())

plt.rc('font', size=16)

def main():
   if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

   pt = os.path.join(ARCHV, CASE, HISTS, H0)
   ds = ux.open_mfdataset(GRIDFN, pt)
   flt = cftime.date2num(ds.time, tunits)
   print(flt)

   print('\nObtaining dp3d...')
   aterm = (ds['hyai'].isel(ilev=slice(1, None)) - ds['hyai'].isel(ilev=slice(None, -1)).data) * P0
   bterm = (ds['hybi'].isel(ilev=slice(1, None)) - ds['hybi'].isel(ilev=slice(None, -1)).data) * ds['PS']
   dp3d = aterm + bterm
   dp3d = dp3d.rename(dict(ilev='lev')).assign_coords(lev=ds['lev'])
   print(dp3d)
   ds = ds.assign(variables=dict(dp3d=dp3d))

   print('\nBasic tests...')
   print('Atmospheric mass from PS')
   print((weighted_temporal_mean(ds['PS']).weighted_mean() * GA / g).values)


   print('Atmospheric mass from 3D integral')
   colintmass = ((ds['U'] * 0 + 1) * dp3d).sum(dim='lev') / g
   print(weighted_temporal_mean(colintmass.weighted_mean()).values * GA)


   for var in VARS:
      print('Working on', var)
      gint, horagg = None, None
      if var in RAW:
         horagg = 'mean'
         colint = (ds[var] * dp3d).sum(dim='lev') / g
         gint = colint.weighted_mean().assign_coords(coords=dict(time=ds['time']))
      else:
         horagg = UDEF[var][-1]
         rawvars = [[ds[var] if type(var) == str else float(var) for var in tup] for tup in UDEF[var][:-1]]
         sumterms = sum([math.prod(rv) for rv in rawvars])
         colint = (sumterms * dp3d).sum(dim='lev') / g
         gint = colint.weighted_mean().assign_coords(coords=dict(time=ds['time']))
         if horagg == 'sum':
            gint *= GA

      with open(os.path.join(OUTDIR, 'globints.txt'), 'a') as f:
         if var == VARS[0]:
            print('\n', pt, file=f)
            print(gint.isel(time=slice(CLIPMO, None)).shape, '\n', file=f)
         print(var, horagg, weighted_temporal_mean(gint.isel(time=slice(CLIPMO, None))), file=f)
         f.close()

      plt.plot(flt, gav.values)
      #plt.legend()
      plt.title(var)
      plt.xlabel(tunits)# [%s]' % str(flt.units))
      plt.ylabel(units)
      plt.savefig(os.path.join(OUTDIR, '%s.png' % var), bbox_inches='tight')
      plt.close()

   print('globint3d.py done.')



#adapted from https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
#weight monthly avgs to correctly compute annual avgs
def weighted_temporal_mean(da):
   """
   weight by days in each month
   """
   # Determine the month length
   month_length = da.time.dt.days_in_month
   #print(month_length)

   num = (da * month_length).sum(dim='time')
   den = month_length.sum(dim='time')
   return num / den

if __name__ == '__main__':
   main()
