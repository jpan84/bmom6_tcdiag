import cftime
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

ciceh = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250206_1degeqm885/ice/hist/*.cice.h.*.nc'
ciceh = '/glade/derecho/scratch/jpan/b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.250130_h80l89/run/cice.r.cat.nc'
grid = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250206_1degeqm885/ice/hist/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250206_1degeqm885.cice.h.0003-02.nc'
tunits = 'common_years since 0000-01-01'
latcoords = 'TLAT'

def main():
   ds = xr.open_mfdataset(ciceh)
   flt = None
   if 'time' in ds.coords:
      flt = cftime.date2num(ds.time, tunits)

   if not latcoords in ds.coords:
      grds = xr.open_dataset(grid)
      ds = ds.rename_dims(dims_dict=dict(record='time', ncat='nc'))
      ds = ds.assign_coords(coords={latcoords: grds[latcoords], 'time': ds.attrs['myear'] + ds.time})
      flt = ds.time

   vol = ds['vicen'].sum(dim='nc')

   nh = vol.isel(nj=slice(366, None))
   sh = vol.isel(nj=slice(0, 100))
   print(nh)
   print(sh)

   nhser = latlon_avg(nh, nh[latcoords])
   shser = latlon_avg(sh, sh[latcoords])
   print(nhser)
   print(shser)

   plt.plot(flt, nhser.values, label='nh')
   plt.plot(flt, shser.values, label='sh')
   plt.xlabel('year')
   plt.ylabel('Mean ice thickness [m]')
   plt.title('Lat: %.1f %.1f' % (nh[latcoords].min().values, nh[latcoords].max().values))
   plt.legend()
   plt.savefig('cice_tser.png')
   plt.show()

def latlon_avg(da, latcoord, latdim='nj', londim='ni'):
   wgts = np.cos(np.deg2rad(latcoord))
   print('weights shape', wgts.shape)
   return ((da * wgts).sum(dim=latdim) / wgts.sum(dim=latdim)).mean(dim=londim)

if __name__ == '__main__':
   main()
