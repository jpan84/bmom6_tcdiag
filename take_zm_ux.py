import os
import sys
import uxarray as ux
import xarray as xr
from dask.diagnostics import ProgressBar

DIRI = '/glade/derecho/scratch/jpan/jpan_tcfields1/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/'
FILI = 'TC_R4_masked_yhoureddy_h1i_0012-0014.nc'
UGRD = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
VARS = ['PS', 'PRECT']
LATS = (-90, 90, 1.)

ds = ux.open_mfdataset(UGRD, os.path.join(DIRI, FILI))
outds = None
for dv in VARS:
   zm = ds[dv].zonal_mean(lats=LATS)
   if outds is None:
      outds = xr.Dataset(data_vars={dv: zm})
   else:
      outds = outds.assign(variables={dv: zm})

outds = outds.assign_attrs(script_from=sys.argv[0])
with ProgressBar():
   outds.to_netcdf(os.path.join(DIRI, FILI + '.zm.nc'))
