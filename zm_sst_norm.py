import sys
import xarray as xr
from dask.diagnostics import ProgressBar

FILI = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/ocn/hist/*.sfc.*.nc'
CASES = ['250415_unseed', '250417_ctrl', '250416_seed1x1']
ALIASES = ['UNSEED', 'CTRL', 'SEED']
FILO = 'zm_sst_ydaymean.nc'

FLDNM = 'tos'
LATNM = 'yh'
LONNM = 'xh'
MINLAT = 7.0
MAXLAT = 16.0

def main():
   print('Opening datasets...')
   dss = [xr.open_mfdataset(FILI % cs) for cs in CASES]

   print('Selecting latitude range...')
   nhsst = [ds[FLDNM].sel(indexers={LATNM: slice(MINLAT, MAXLAT)}) for ds in dss]
   shsst = [ds[FLDNM].sel(indexers={LATNM: slice(-MAXLAT, -MINLAT)}) for ds in dss]
   selsst = [xr.concat([nhsst[ii], shsst[ii]]).expand_dims(case=[ALIASES[ii]]) for ii in range(len(nhsst))]

   print('Taking zonal and day-of-year means...')
   zmnorm = xr.concat([ss.mean(dim=LONNM).groupby('time.dayofyear').mean() for ss in selsst], dim='case')

   print('Saving .nc...')
   with ProgressBar():
      xr.Dataset(data_vars={FLDNM: zmnorm}).to_netcdf(FILO)

   print(sys.argv[0], 'done.')

if __name__ == '__main__':
   main()
