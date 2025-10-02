import uxarray as ux
import xarray as xr
import numpy as np
import os

DIRI = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/'
VTFIL = 'VT500_yhoureddy_0012-0014.nc'
VUFIL = 'VU500_yhoureddy_0012-0014.nc'
CAMGRID = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

LATBND = 30
LATS = (-LATBND, LATBND, 2)

def main():
   EMFds = ux.open_mfdataset(CAMGRID, os.path.join(DIRI, VUFIL))
   EHFds = ux.open_mfdataset(CAMGRID, os.path.join(DIRI, VTFIL))
   szn = EMFds['time'].dt.season

   EMF500_zm = EMFds['V500'].zonal_mean(lat=LATS)
   EHF500_zm = EHFds['V500'].zonal_mean(lat=LATS)
   EMF500_zm, EHF500_zm = EMF500_zm.assign_coords(time=EMFds.time), EHF500_zm.assign_coords(time=EHFds.time)
   print(EMF500_zm.sel(time=szn.isin(['DJF', 'MAM']), latitudes=slice(0, EMF500_zm.latitudes[-1])))

   zmds = xr.Dataset(data_vars=dict(EHF500=EHF500_zm.to_xarray(), EMF500=EMF500_zm.to_xarray()))
   nhwarmds = zmds.sel(time=szn.isin(['JJA', 'SON']), latitudes=slice(-5, LATBND))
   shwarmds = zmds.sel(time=szn.isin(['DJF', 'MAM']), latitudes=slice(-LATBND, 5))
   #shmir = shwarmds.isel(latitudes=slice(None, None, -1)).assign_coords(latitudes=nhwarmds.latitudes)

   ###too slow to save shmir
   #outds = xr.concat([nhwarmds, shmir], dim='time')
   #outds = outds.map(lambda da: da.assign_attrs(dict(zonal_mean='True')))
   #outds.to_netcdf(os.path.join(DIRI, 'EHF_EMF_500_warmhemi_0012-0014.nc'))

   saveds = lambda ds, hemi: ds.map(lambda da: da.assign_attrs(dict(zonal_mean='True')))\
                     .to_netcdf(os.path.join(DIRI, 'EHF_EMF_500_warm%s_0012-0014.nc' % hemi))
   saveds(nhwarmds, 'NH')
   saveds(shwarmds, 'SH')

if __name__ == '__main__':
   main()
