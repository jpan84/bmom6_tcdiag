import xarray as xr
import numpy as np
import os
import sys

DIR = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_regrid_0.25x0.25_onpres'

ds = xr.open_mfdataset(os.path.join(DIR, '*.h0a.*.nc'))
print(ds.time.dt.days_in_month.data)
dpm = ds.time.dt.days_in_month.data
wgts = dpm / dpm.sum()

for dv in ds.data_vars:
   if ds[dv].dims[0] == 'time' and np.issubdtype(ds[dv].dtype, np.number):
      print('Averaging', dv)
      ltm = np.einsum('i,i...->...', wgts, ds[dv]) #output of to_dataarray() has 'variable' as new 1st dimension. CESM output has time as 1st dimension.
      ltm = xr.DataArray(ltm, coords={dm: ds[dm] for dm in ds[dv].dims[1:]}, attrs=ds[dv].attrs)
      ds = ds.assign(variables={dv: ltm})
   else:
      print('Skipping', dv)

ds.attrs['history'] = sys.argv[0] + ';\n' + ds.attrs['history']
ds.attrs['ltm_start_date'] = str(ds.time.values[0])
ds.attrs['ltm_end_date'] = str(ds.time.values[-1])
ds.to_netcdf(os.path.join(DIR, 'ltm.nc'))
