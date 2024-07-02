#Joshua Pan Feb 2024
#compute eofs of monthly SST (tos)

import sys
import os
import numpy as np
import xarray as xr
from eofs.xarray import Eof
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs

LATLAB = np.array([-90., -60., -30., 0., 30., 60., 90.])
latrng = (-90, 90)
OUTDIR = 'sstvar/'
hists = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.f09_sx0.66av1.aqua.production.0530ywbranch/ocn/hist/*mom6.hm_*.nc'
DV = 'tos'

def main():
   print('Opening datasets...')
   dss = xr.open_mfdataset(hists, chunks={'time': 12})

   print('Selecting...')
   sst = ds[DV].sel(yh=slice(*latrng))
   sst = sst.copy(data=sst.values) #copy array so that eofs doesn't have to deal with dask

   print('Solving EOFs...')
   coslat = np.cos(np.deg2rad(sst.yh.values)).clip(0., 1.)
   wgts = np.sqrt(coslat)
   solver = Eof(sst, weights=wgts)

   print('Extracting EOFs and PCs...')
   eofs = solver.eofsAsCovariance(pcscaling=1, neofs=6) #pcs scaled to unit variance
   pcs = solver.pcs(pcscaling=1, npcs=6)
   varfracs = solver.varianceFraction(neigs=10)

   print('Saving EOFs...')
   outds = xr.Dataset(data_vars=dict(eofs=eofs, pcs=pcs, expvar=varfracs))
   outds.to_netcdf(os.path.join(OUTDIR, '%s_ssteofs.nc' % hists.split('/')[-4]))

   '''
   print('Aligning signs...')
   for ii in range(len(eofs)):
      for nn in eofs[ii].mode:
         ff = eofs[ii].sel(mode=nn)
         if ff.max() < -ff.min():
            eofs[ii].loc[dict(mode=nn)] = -eofs[ii].sel(mode=nn)
            pcs[ii].loc[dict(mode=nn)] = -pcs[ii].sel(mode=nn)
   '''

   print('Plotting EOFs and PCs...')
   plt.rc('font', size=20)
   for nn in eofs.mode:
      ax = plt.axes(projection=ccrs.PlateCarree())
      cset = ax.contourf(eofs.xq, eofs.yh, eofs.sel(mode=nn).values, cmap='RdBu_r')
      ax.coastlines()
      plt.colorbar(cset)
      plt.savefig(os.path.join(OUTDIR, 'eof%d.png' % nn), bbox_inches='tight')
      plt.close()

      plt.plot(pcs.time - pcs.time[0], pcs.sel(mode=nn).values)
      plt.legend()
      plt.xlabel('Month')
      plt.ylabel('pc%d' % nn)
      plt.savefig(os.path.join(OUTDIR, 'pc%d.png' % nn), bbox_inches='tight')
      plt.close()

   print('Plotting explained variance...')
   plt.plot(varfracs.mode, varfracs.values)
   plt.xlabel('Mode')
   plt.ylabel('Explained variance')
   plt.legend()
   plt.savefig(os.path.join(OUTDIR, 'expvar.png'), bbox_inches='tight')
   plt.close()

   print('%s done.' % sys.argv[0])

if __name__ == '__main__':
   main()
