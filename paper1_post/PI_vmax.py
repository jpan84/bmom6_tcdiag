import sys
import os
import uxarray as ux
import xarray as xr
import numpy as np
from tcpyPI import pi as tcpi
from paths import CAMGR, ARCHRT, ALIA
from consts import TCK, OM
from sznl_funcs import monthly2sznl, stack_hemi_sznl
import matplotlib.pyplot as plt

from dask.distributed import Client

TSTFIL = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist_onpres_gnupar/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.cam.h0a.0011-08.onpres.nc'
dates = '000[8-9]'
hist_onp = 'atm/hist_onpres_gnupar/*.h0a.%s*.onpres.nc' % dates

FILO = os.path.join(ARCHRT[0], 'atm/', 'vmax_PI_ne120_%s.nc' % dates)

pintp = lambda x, p: x.interp(plev=p).drop('plev')

def main():
   # This automatically spins up a local cluster and hooks xarray/uxarray into it
   client = Client() 
   print(f"Dashboard available at: {client.dashboard_link}")

   dss = [ux.open_mfdataset(CAMGR, os.path.join(ar, hist_onp)).expand_dims(case=[ALIA[ii]]) for ii, ar in enumerate(ARCHRT)]
   for ds in dss:
      ds.uxgrid = dss[0].uxgrid
   ds = ux.concat(dss, 'case')
   ds = ds.sortby('plev', ascending=False)
   ds = ds.chunk({'time': 1, 'plev': -1, 'n_face': -1}) #correct chunking provided the biggest speedup (several hours down to 0.5 hr for 1 year)
   print(ds)

   print('Task graphs for PI args...')
   SSTC = ds['TS'] - TCK
   MSL = ds['PSL'] / 100
   phPa = ds['plev'] / 100
   TC = ds['T'] - TCK
   rwv = ds['Q'] / (1 - ds['Q']) * 1000

   #PI vmax
   print('PI task graph...')
   vmax, pmin, ifl, _, _ = pi_xr_ux(SSTC, MSL, phPa, TC, rwv)
   vmax = ux.UxDataArray(vmax, uxgrid=ds.uxgrid)
   pmin = ux.UxDataArray(pmin, uxgrid=ds.uxgrid)

   print('Saving to netcdf...')
   # 1. Convert vmax to a clean, standard numpy or dask array stripped of parent metadata
   vmax_clean = xr.DataArray(
       vmax.data, 
       dims=['case', 'time', 'n_face'],
       coords={
         'case': ds.case.values,
         'time': ds.time.values
       }
   )
   pmin_clean = xr.DataArray(
       pmin.data, 
       dims=['case', 'time', 'n_face'],
       coords={
         'case': ds.case.values,
         'time': ds.time.values
       }
   )
   xr.Dataset(data_vars=dict(vmax_PI=vmax_clean, pmin_PI=pmin_clean)).to_netcdf(FILO)
   print('Saved.')
   exit()

def main_plot():
   ds = ux.open_mfdataset(CAMGR, FILO)
   print(ds['vmax_PI'].values)

   print('Zonal mean task graph...')
   vmax_zm = ds['vmax_PI'].zonal_mean(lat=(-40, 40, 2)).squeeze()
   vmax_zm = xr.DataArray(vmax_zm).assign_coords(coords=ds['vmax_PI'].coords)

   print('Time agg task graph...')
   vmax_mm = vmax_zm.groupby('time.month').mean()
   vmax_sm = stack_hemi_sznl(monthly2sznl(vmax_mm)).mean(dim='season')

   print(vmax_sm.values)

   print('Computing and plotting...')
   [plt.plot(vmax_zm.latitudes, vmax_sm.sel(case=cs), label=str(cs.values)) for ii, cs in enumerate(vmax_zm['case'])]
   plt.legend()
   plt.show()


def pi_xr_ux(*piargs):
   return xr.apply_ufunc(tcpi, *piargs,\
            input_core_dims=[[], [], ['plev'], ['plev'], ['plev']],\
            output_core_dims=[[] for _ in range(len(piargs))],\
            vectorize=True, dask='parallelized',\
            output_dtypes=[float, float, int, float, float])

if __name__ == '__main__':
   if sys.argv[1] == 'compute':
      main()
   if sys.argv[1] == 'plot':
      main_plot()
