import os
import sys
import xarray as xr
from dask.diagnostics import ProgressBar

#FILI = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/ocn/hist/*.sfc.*.nc'
#CASES = ['250415_unseed', '250417_ctrl', '250416_seed1x1']
pths = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250702_unseed2hPa6m/ocn/hist/', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/ocn/hist', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/ocn/hist', '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/ocn/hist', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/ocn/hist/']
ALIASES = ['UNSEED2', 'UNSEED', 'CTRL', 'MSEED', 'SEED']
FILO = 'zm_sst.nc'
FILONORM = 'sstmaxlat_ydaymean.nc'

FLDNM = 'tos'
LATNM = 'yh'
LONNM = 'xh'
MINLAT = 7.0
MAXLAT = 16.0

def main():
   print('Opening datasets...')
   dss = [xr.open_mfdataset(os.path.join(pt, '*.sfc.*.nc')) for pt in pths]

   print('Selecting latitude range...')
   #nhsst = [ds[FLDNM].sel(indexers={LATNM: slice(MINLAT, MAXLAT)}) for ds in dss]
   #shsst = [ds[FLDNM].sel(indexers={LATNM: slice(-MAXLAT, -MINLAT)}) for ds in dss]
   #selsst = [xr.concat([nhsst[ii], shsst[ii]], dim=LATNM).expand_dims(case=[ALIASES[ii]]) for ii in range(len(nhsst))]
   selsst = [ds[FLDNM].sel(indexers={LATNM: slice(-MAXLAT, MAXLAT)}) for ds in dss]

   print('Taking zonal means...')
   zmsst = xr.concat([ss.mean(dim=LONNM) for ss in selsst], dim='case')

   print('Taking zonal and day-of-year means...')
   #zmnorm = zmsst.groupby('time.dayofyear').mean()
   #zmnorm = xr.concat([ss.mean(dim=LONNM).groupby('time.dayofyear').mean() for ss in selsst], dim='case')
   zmmaxlat = zmsst.idxmax(LATNM)
   latnorm = zmmaxlat.groupby('time.dayofyear').mean()

   print('Saving .nc files...')
   with ProgressBar():
      xr.Dataset(data_vars={FLDNM: zmsst}).to_netcdf(FILO)
      #xr.Dataset(data_vars={FLDNM: zmnorm}).to_netcdf(FILONORM)
      xr.Dataset(data_vars={LATNM: latnorm}).to_netcdf(FILONORM)

   print(sys.argv[0], 'done.')

if __name__ == '__main__':
   main()
