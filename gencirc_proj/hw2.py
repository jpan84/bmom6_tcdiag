import stereoplot
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

FILI = '/glade/work/jpan/597_gencirc_gfdl/model_results/postprocessed_251003_jpan_DJF_zonsym_midlat.nc'
QFIL = '/glade/work/jpan/597_gencirc_gfdl/ModelCopy/input/midlatitude_forcing_from_HK81.nc'
DIRO = 'hw2_figs/'

def main():
   ds = xr.open_dataset(FILI, decode_times=False)
   qds = xr.open_dataset(QFIL)

   print(qds['heating'])
   #print(qds['heating'].sel(level=0.5, method='nearest').max())

   plt_psi = ds['psi_anom'].isel(time=slice(14, None)).mean(dim='time')
   plt_psi = plt_psi.sel(level=850)
   #plt_psi = (plt_psi.sel(level=250) - plt_psi.sel(level=850)) / 2
   psilevs = 1e6 * np.arange(1, 21, 1)
   psilevs = np.concatenate((-psilevs[::-1], psilevs))
   print(ds.time.isel(time=slice(14, None)))

   ax = stereoplot.NPstereoaxes(0)
   ax.contour(qds['longitude'], qds['latitude'], qds['heating'].sel(level=0.5, method='nearest') * 86400, transform=ccrs.PlateCarree(), levels=[1], colors='red')
   ax.contour(plt_psi['longitude'], plt_psi['latitude'], plt_psi, transform=ccrs.PlateCarree(), levels=psilevs, colors='black')
   ax.set_title('(b)', loc='left')

   plt.show()
   plt.close()

   plot_u = ds['u_anom'].isel(time=slice(14, None)).mean(dim=['time', 'longitude'])
   plt.contourf(plot_u['latitude'], plot_u['level'], plot_u, cmap='bwr', levels=np.arange(-2, 2.1, .25))
   plt.xlim(-10, 90)
   plt.xlabel('lat')
   plt.ylim(1000, 100)
   plt.yscale('log')
   plt.ylabel('pressure [hPa]')
   plt.colorbar()
   plt.show()

if __name__ == '__main__':
   main()
