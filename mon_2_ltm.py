import xarray as xr
import numpy as np
import os
import sys

DIR = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_regrid_0.25x0.25_onpres'

ds = xr.open_mfdataset(os.path.join(DIR, '*.h0a.*.nc'))
#print(ds.time.dt.days_in_month.data)
dpm = ds.time.dt.days_in_month.data
wgts = dpm / dpm.sum()

fullyrs = wgts.size // 12

for ii in range(fullyrs):
   tslc = slice(12 * ii, 12 * (ii + 1))
   subds = ds.isel(time=tslc)
   subwgts = wgts[tslc]

   for dv in subds.data_vars:
      if subds[dv].dims[0] == 'time' and np.issubdtype(subds[dv].dtype, np.number):
         print('Averaging', dv)
         ltm = np.einsum('i,i...->...', subwgts, subds[dv]) #output of to_dataarray() has 'variable' as new 1st dimension. CESM output has time as 1st dimension.
         ltm = xr.DataArray(ltm, coords={dm: subds[dm] for dm in subds[dv].dims[1:]}, attrs=subds[dv].attrs)
         ltm = ltm.assign_coords(coords=dict(time=subds.time.values[0]))
         subds = subds.assign(variables={dv: ltm})
      else:
         print('Skipping', dv)
   
   subds.attrs['history'] = sys.argv[0] + ';\n' + subds.attrs['history']
   subds.attrs['ltm_start_date'] = str(subds.time.values[0])
   subds.attrs['ltm_end_date'] = str(subds.time.values[-1])
   subds.to_netcdf(os.path.join(DIR, 'ltm_%d.nc' % ii))
