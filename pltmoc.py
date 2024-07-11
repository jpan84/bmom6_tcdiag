#plot merid overturning resid streamf based on https://mom6-analysiscookbook.readthedocs.io/en/latest/notebooks/Overturning.html
#import warnings
#warnings.filterwarnings('ignore')

import os
import numpy as np
import xarray as xr # requires >= 0.15.1
import pandas as pd
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt
from xoverturning import calcmoc

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0530ywbranch/'
HISTS = 'ocn/hist/'
H0 = r'*mom6.hm_[0-9]*.nc'
diro = 'sstvar/'

def main():
   print('Opening history files...')
   ds = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H0)).rename(dict(zl='z_l', zi='z_i'))
   dsgrid = xr.open_dataset('/glade/derecho/scratch/jpan/b.e23.BMOM.f09_sx0.66av1.aqua.production.0515ctrl/run/INPUT/sx0.66av1_grid.nc')
   geo = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, '*static*.nc'))

   print(ds.dims)
   print(dsgrid.dims)

   moc = calcmoc(xr.merge([ds, geo]))#, dsgrid=dsgrid)
   moc = moc.rename('psieul')
   print(moc)
   #outds = xr.Dataset(data_vars=dict(psieul=moc))
   #outds.to_netcdf(os.path.join(diro, '%s_moc.nc' % CASE[:-1].split('.')[-1]))
   #print(moc.time.values)
   print('Saving moc to nc...') 
   towrt = moc.to_netcdf(path=os.path.join(diro, '%s_moc.nc' % CASE[:-1].split('.')[-1]), compute=False)
   with ProgressBar():
      towrt.compute()

   moc_mean = None
   with ProgressBar():
      moc_mean = moc.isel(time=slice(0,60)).mean('time').load()

   plt.contourf(moc.yq, moc.z_i, moc_mean.values, levels=np.arange(-50, 51, 5), cmap='RdBu_r', extend='both')
   plt.gca().invert_yaxis()
   plt.xlabel('lat')
   plt.ylabel('z [m]')
   plt.colorbar()
   plt.savefig('%s_psieul.png' % CASE[:-1].split('.')[-1])
   plt.show()

if __name__ == '__main__':
   main()
