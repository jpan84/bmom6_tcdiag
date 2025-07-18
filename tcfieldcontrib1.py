#compute the fraction or component of an h1i field that is accounted for by TCs
#Use raw h1i output and NodeFileFilter'ed h1i output, and yhour eddy fields

import sys
import os
import operator
import uxarray as ux
import xarray as xr
import numpy as np
from sznl_funcs import monthly2sznl, stack_hemi_sznl
import pltsettings
import matplotlib.pyplot as plt


ARCHV = '/glade/campaign/univ/upsu0032/jpan_tcfields/'
ALIAS = '250417_ctrl'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s' % ALIAS
DIRO = './tcfields2mps_%s/' % ALIAS
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

RAWALL = 'hist_0010_h1i/cat_h1i.nc'
RAWTCS = 'TC_R2_masked/*.h1i.*.nc'
EDDALL = 'yhoureddy/yhoureddy_h1i.nc'
EDDTCS = 'yhoureddy_TC_R2_masked/yhoureddy_h1i.nc'

zmlats = (-90, 90.1, 0.5)
LATLAB = np.arange(-90, 91, 30)
g = 9.81

ZMLATS = np.arange(*zmlats)
SINLAT = np.sin(np.deg2rad(ZMLATS))
LSTYS = ['solid', 'dashed']
LCLRS = ['blue', 'orange']

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)
   pltsettings.set1()

   print('Opening datasets...')
   rawallds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, RAWALL))
   rawtcsds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, RAWTCS))
   eddallds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, EDDALL))
   eddtcsds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, EDDTCS))

   outallds, outtcsds = compute_unsigned_flds(rawallds, rawtcsds, ['PRECT', 'TAUX'])

   print('Saving .nc outputs...')
   outallds.to_netcdf(os.path.join(DIRO, 'means_all.nc'))
   outtcsds.to_netcdf(os.path.join(DIRO, 'means_tcs.nc'))

   print(sys.argv[0], 'done.')

#Compute seasonal and zonal mean after all physical quantity computation/masking is complete
def all_and_TC_to_sznlzm(allda, tcsda):
   monzm = [da.groupby('time.month').mean().zonal_mean(lat=ZMLATS) for da in [allda, tcsda]]
   return [stack_hemi_sznl(monthly2sznl(da)) for da in monzm]

def compute_unsigned_flds(allds, tcsds, flds):
   retall, rettcs = None, None
   for fl in flds:
      print('Processing %s unsigned...' % fl)
      sznzm = all_and_TC_to_sznlzm(allds[fl], tcsds[fl])
      if retall is None:
         retall = xr.Dataset(data_vars={fl: sznzm[0].to_xarray()})
         rettcs = xr.Dataset(data_vars={fl: sznzm[1].to_xarray()})
      else:
         retall = retall.assign(variables={fl: sznzm[0].to_xarray()})
         rettcs = rettcs.assign(variables={fl: sznzm[1].to_xarray()})
   return retall, rettcs


if __name__ == '__main__':
   main()
