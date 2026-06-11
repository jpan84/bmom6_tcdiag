import uxarray as ux
import xarray as xr
from tcpyPI import pi as tcpi
from paths import CAMGR
from consts import TCK
import matplotlib.pyplot as plt

TSTFIL = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist_onpres_gnupar/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.cam.h0a.0011-08.onpres.nc'


def main():
   ds = ux.open_mfdataset(CAMGR, TSTFIL)
   ds = ds.sortby('plev', ascending=False)
   ds = ds.chunk({'plev': -1})
   print(ds)

   SSTC = ds['TS'] - TCK
   MSL = ds['PSL'] / 100
   phPa = ds['plev'] / 100
   TC = ds['T'] - TCK
   rwv = ds['Q'] / (1 - ds['Q']) * 1000

   vmax, _, ifl, _, _ = pi_xr_ux(SSTC, MSL, phPa, TC, rwv)
   vmax = ux.UxDataArray(vmax, uxgrid=ds.uxgrid)
   vmax_zm = vmax.zonal_mean(lat=(-50, 50, 2)).squeeze()
   print('IFL', ifl.values)
   print(vmax_zm.values)

   plt.plot(vmax_zm.latitudes, vmax_zm)
   plt.show()

def pi_xr_ux(*piargs):
   return xr.apply_ufunc(tcpi, *piargs,\
            input_core_dims=[[], [], ['plev'], ['plev'], ['plev']],\
            output_core_dims=[[] for _ in range(len(piargs))],\
            vectorize=True, dask='parallelized',\
            output_dtypes=[float, float, int, float, float])

if __name__ == '__main__':
   main()
