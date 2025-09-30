#Ex: python3 -u hyb_delta_pressure.py /glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist dp_monthly.nc
import sys
import os
import xarray as xr
from dask.diagnostics import ProgressBar

DIRI = sys.argv[1]
tape = '*h0a*.nc'
outnm = sys.argv[2]

P0 = 1e5 #Pa

def main():
   print('Opening dataset...')
   ds = xr.open_mfdataset(os.path.join(DIRI, tape))

   print('\"Computing\" dp...')
   aterm = (ds['hyai'].isel(ilev=slice(1, None)) - ds['hyai'].isel(ilev=slice(None, -1)).data) * P0
   bterm = (ds['hybi'].isel(ilev=slice(1, None)) - ds['hybi'].isel(ilev=slice(None, -1)).data) * ds['PS']
   dp3d = aterm + bterm
   dp3d = dp3d.rename(dict(ilev='lev')).assign_coords(lev=ds['lev'])

   print('Saving dp...')
   with ProgressBar():
      xr.Dataset(data_vars=dict(dp3d=dp3d)).to_netcdf(os.path.join(DIRI, outnm))

if __name__ == '__main__':
   main()
