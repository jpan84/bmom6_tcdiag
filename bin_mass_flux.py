import numpy as np
from scipy.stats import binned_statistic_2d
#import uxarray as ux
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

g = 9.81
cp = 1004
lv = 2500840
a = 6.371e6

testfile = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_0010_h1i/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.0010-*-01-00000.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

#TODO: allow flexible latitude bins?
#TODO: allow diffing cases
#TODO: split seasons
def main():
   ds = xr.open_mfdataset(testfile)
   print(ds)
   #ds = ds.isel(ncol=(np.abs(ds['lat']) < 35).all(dim='time').compute())
   nh, sh = sel_unstruct_tropics(ds)
   print(nh)
   print(sh)
   ds = xr.concat((nh, sh), 'ncol')

   mse850 = g * ds['Z850'] + cp * ds['T850'] + lv * ds['Q850']
   plt.hist(mse850.data)
   plt.savefig('mse850_hist.png')
   plt.close()

   area = ds['area'] * a**2
   umf500 = (ds['OMEGA500'] < 0) * (-ds['OMEGA500'] * area / g)

   msebins = np.arange(3e5, 3.61e5, 5e3)
   capebins = np.arange(0, 501, 25)

   umf_b_2d, mse_edges, cape_edges, _ = bin_umf_2d(umf500, mse850, ds['CAPE'], msebins, capebins, ds['time'].size)
   umf_b_2d = np.nan_to_num(umf_b_2d, nan=0.0)
   print(mse_edges, cape_edges, umf_b_2d.shape)
   plt.pcolormesh(mse_edges, cape_edges, umf_b_2d.T, shading='flat')#, norm=colors.LogNorm())
   plt.xlabel('MSE 850 [J kg$^{-1}$]')
   plt.ylabel('CAPE [J kg$^{-1}$]')
   plt.title('500 hPa UMF [kg s$^{-1}$]')
   plt.colorbar()
   plt.savefig('umf_binned_2d.png', bbox_inches='tight')
   plt.close()

def sel_unstruct_tropics(ds, outerlat=30., hemi='warm'):
   szn = ds['time'].dt.season
   if hemi == 'warm':
      nh = ds.sel(time=szn.isin(['JJA', 'SON'])).isel(ncol=((ds['lat'] > 0) & (ds['lat'] < outerlat)).all(dim='time').compute())
      sh = ds.sel(time=szn.isin(['DJF', 'MAM'])).isel(ncol=((ds['lat'] < 0) & (ds['lat'] > -outerlat)).all(dim='time').compute())
      return nh, sh
   if hemi == 'cool':
      nh = ds.sel(time=szn.isin(['DJF', 'MAM'])).isel(ncol=((ds['lat'] > 0) & (ds['lat'] < outerlat)).all(dim='time').compute())
      sh = ds.sel(time=szn.isin(['JJA', 'SON'])).isel(ncol=((ds['lat'] < 0) & (ds['lat'] > -outerlat)).all(dim='time').compute())
      return nh, sh


def bin_umf_2d(umf, varx, vary, binsx, binsy, ntime):
   umfsum, xedges, yedges, binnum = binned_statistic_2d(x=varx.data.ravel(), y=vary.data.ravel(),\
                 values=umf.data.ravel(), statistic='sum', bins=[binsx, binsy])
   return umfsum / ntime, xedges, yedges, binnum

if __name__ == '__main__':
   main()
