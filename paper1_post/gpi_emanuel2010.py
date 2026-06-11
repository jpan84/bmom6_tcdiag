import uxarray as ux
import xarray as xr
from tcpyPI import pi as tcpi
from paths import CAMGR
from consts import TCK

TSTFIL = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist_onpres_gnupar/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.cam.h0a.0011-08.onpres.nc'


def main():
   ds = ux.open_mfdataset(CAMGR, TSTFIL)
   ds = ds.chunk({"plev": -1})
   print(ds)

   SSTC = ds['TS'] - TCK
   MSL = ds['PSL'] / 100
   phPa = ds['plev'] / 100
   TC = ds['T'] - TCK
   rwv = ds['Q'] / (1 - ds['Q']) * 1000

   print(pi_xr_ux(SSTC, MSL, phPa, TC, rwv))

def pi_xr_ux(*piargs):
   return xr.apply_ufunc(tcpi, *piargs,\
            input_core_dims=[[], [], ['plev'], ['plev'], ['plev']],\
            output_core_dims=[[] for _ in range(len(piargs))],\
            vectorize=True, dask='parallelized',\
            output_dtypes=[float, float, int, float, float])

if __name__ == '__main__':
   main()
