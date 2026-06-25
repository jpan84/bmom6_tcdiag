import os
import uxarray as ux
import xarray as xr
import numpy as np
from tcpyPI import pi as tcpi
from paths import CAMGR, ARCHRT, ALIA
from consts import TCK, OM
import matplotlib.pyplot as plt

TSTFIL = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist_onpres_gnupar/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.cam.h0a.0011-08.onpres.nc'
hist_onp = 'atm/hist_onpres_gnupar/*.onpres.nc'

pintp = lambda x, p: x.interp(plev=p).drop('plev')

def main():
   dss = [ux.open_mfdataset(CAMGR, os.path.join(ar, hist_onp)).expand_dims(case=[ALIA[ii]]) for ii, ar in enumerate(ARCHRT)]
   ds = ux.concat(dss, 'case')
   ds = ds.sortby('plev', ascending=False)
   ds = ds.chunk({'plev': -1})
   print(ds)

   print('Task graphs for PI args...')
   SSTC = ds['TS'] - TCK
   MSL = ds['PSL'] / 100
   phPa = ds['plev'] / 100
   TC = ds['T'] - TCK
   rwv = ds['Q'] / (1 - ds['Q']) * 1000

   #PI vmax
   print('PI task graph...')
   vmax, _, ifl, _, _ = pi_xr_ux(SSTC, MSL, phPa, TC, rwv)
   vmax = ux.UxDataArray(vmax, uxgrid=ds.uxgrid)

   print('Saving to netcdf...')
   xr.Dataset(data_vars=dict(vmax_PI=vmax)).to_netcdf(os.path.join(ARCHRT[0], 'atm/', 'vmax_PI_ne120.nc'))
   exit()

   print('Zonal mean task graph...')
   vmax_zm = vmax.zonal_mean(lat=(-40, 40, 2)).squeeze()
   #print('IFL', ifl.values)
   #print(vmax_zm.values)

   #plt.plot(vmax_zm.latitudes, vmax_zm)
   #plt.show()

   print('Computing and plotting...')
   plt.plot(gpi_zm.latitudes, gpi_zm)
   plt.show()


def pi_xr_ux(*piargs):
   return xr.apply_ufunc(tcpi, *piargs,\
            input_core_dims=[[], [], ['plev'], ['plev'], ['plev']],\
            output_core_dims=[[] for _ in range(len(piargs))],\
            vectorize=True, dask='parallelized',\
            output_dtypes=[float, float, int, float, float])

if __name__ == '__main__':
   main()
