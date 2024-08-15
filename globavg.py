#Joshua Pan Oct 2023
#compute globally area-weighted averaged quantities from unstructured grid

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.ne120pg3_sx0.66av1.aqua.production.0812/'
HISTS = 'atm/hist/'
H0 = r'*h0a.[0-9]*.nc'
GRIDDIR = '/glade/p/cesmdata/inputdata/share/scripgrids'
GRIDFN = 'ne120pg3_scrip_170628.nc'
OUTDIR = './globavtraces_0812'

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

vrs = ['FSNT', 'FLNT', 'PRECC', 'PRECL', 'Q', 'QRS', 'QRL', 'PS', 'TMQ', 'QFLX', 'PSDRY', 'PRECT']

def main():
   pt = os.path.join(ARCHV, CASE, HISTS, H0)
   ds = ux.open_mfdataset(os.path.join(GRIDDIR, GRIDFN), pt)

   '''
   #pltsettings.set1()
   args = [(TOPDIR, DIR1, FN), (TOPDIR, DIR2, FN), (TOPDIR3, DIR3, FN)]
   aliases = ['seed', 'unseed', 'ctrl']
   paths = [sorted(glob.glob(os.path.join(*a))) for a in args]
   dss = [ux.open_mfdataset(os.path.join(GRIDDIR, GRIDFN), ps) for ps in paths]
   '''

   for var in vrs:
      gav = globav(ds, var)
      with open(os.path.join(OUTDIR, 'globavgs.txt'), 'a') as f:
         if var == vrs[0]:
            print(args, '\n', file=f)
         print(var, tav(gav), file=f)
         f.close()
      plt.plot(gav.time, gav.values)
      plt.legend()
      plt.title(var)
      plt.xlabel('Time [%s]' % str(gav.time.units))
      plt.ylabel(ds[var].units)
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

if __name__ == '__main__':
   main()
