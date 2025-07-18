#compute the fraction or component of an h1i field that is accounted for by TCs
#Use raw h1i output and NodeFileFilter'ed h1i output, and yhour eddy fields

import sys
import os
import math
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

g = 9.81
unsigned_vars = ['PRECT']
signed_vars = dict(SHFLX=(1, 'SHFLX', False), LHFLX=(1, 'LHFLX', False), TAUX=(1, 'TAUX', False),\
                   TAUY=(1, 'TAUY', True), UMF500=(-1/g, 'OMEGA500', False), UMF850=(-1/g, 'OMEGA850', False),
                   vEHF850=('V850', 'T850', True), vEMF200=('V200', 'U200', True), wEHF500=('OMEGA500', 'U500', False)) #template: (...take the product of these scalars/data_vars..., antisym?)

zmlats = (-90, 90.1, 0.5)
LATLAB = np.arange(-90, 91, 30)

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

   outallds, outtcsds = compute_unsigned_flds(rawallds, rawtcsds, unsigned_vars)

   sgnallds, sgntcsds = compute_signed_flds(rawallds, rawtcsds, signed_vars)
   outallds, outtcsds = xr.merge([outallds, sgnallds]), xr.merge([outtcsds, sgntcsds])

   print('Saving .nc outputs...')
   outallds.to_netcdf(os.path.join(DIRO, 'means_all.nc'))
   outtcsds.to_netcdf(os.path.join(DIRO, 'means_tcs.nc'))

   print(sys.argv[0], 'done.')

#compute the product of terms (scalar or var name) in template
def template_prod(ds, templ):
   terms = [ds[var] if type(var) == str else float(var) for var in templ]
   return math.prod(terms)

#Compute seasonal and zonal mean after all physical quantity computation/masking is complete
def all_and_TC_to_sznlzm(allda, tcsda, antisym=False):
   monzm = [da.groupby('time.month').mean().zonal_mean(lat=ZMLATS) for da in [allda, tcsda]]
   return [stack_hemi_sznl(monthly2sznl(da), antisym=antisym) for da in monzm]

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

def compute_signed_flds(allds, tcsds, flds_dict):
   retall, rettcs = None, None
   for fl in flds:
      print('Processing %s signed...' % fl)
      allfld = template_prod(allds, signed_vars[fl][:-1])
      tcsfld = template_prod(tcsds, signed_vars[fl][:-1])

      for sgn in [('_pos', operator.gt), ('_neg', operator.lt)]:
         allsgn = (sgn[1](allfld, 0)) * allfld
         tcssgn = (sgn[1](tcsfld, 0)) * tcsfld
         sznzm = all_and_TC_to_sznlzm(allsgn, tcssgn, antisym=signed_vars[fl][-1])
         if retall is None:
            retall = xr.Dataset(data_vars={fl + sgn[0]: sznzm[0].to_xarray()})
            rettcs = xr.Dataset(data_vars={fl + sgn[0]: sznzm[1].to_xarray()})
         else:
            retall = retall.assign(variables={fl + sgn[0]: sznzm[0].to_xarray()})
            rettcs = rettcs.assign(variables={fl + sgn[0]: sznzm[1].to_xarray()})
   return retall, rettcs


if __name__ == '__main__':
   main()
