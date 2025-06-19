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

testfile = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_0010_h1i/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.0010-07-01*.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

def main():
   #ds = ux.open_mfdataset(camgrid, testfile)
   ds = xr.open_mfdataset(testfile)
   ds = ds.isel(ncol=(np.abs(ds['lat']) < 35).all(dim='time').compute())
   #print(np.where(np.abs(ds['lat']) < 35)[0].shape)
   #print(ds['area'].min().values, ds['area'].max().values) 
   print(ds)

   mse850 = g * ds['Z850'] + cp * ds['T850'] + lv * ds['Q850']
   plt.hist(mse850.data)
   plt.savefig('mse850_hist.png')
   plt.close()

   area = ds['area'] * a**2
   umf500 = (ds['OMEGA500'] < 0) * (-ds['OMEGA500'] * area / g)

   msebins = np.arange(3e5, 3.61e5, 5e3)
   ###msebinidx = np.digitize(np.ravel(mse850.data), msebins)

   ####serial way
   ######umf_binned = np.empty_like(msebins)
   ######for ii in range(umf_binned.size):
   ######   inbin = np.where(msebinidx == ii)[0]
   ######   umf_binned[ii] = umf500.data.ravel()[inbin].sum()

   ####vectorized
   ###umf_binned = np.bincount(msebinidx, weights=umf500.data.ravel(), minlength=msebins.size) / ds['time'].size

   #umf_b_h850 = bin_umf_by_var1(umf500, mse850, msebins, ds['time'].size)
   #plt.plot(msebins, umf_b_h850)
   #plt.savefig('umf_binned_mse850.png')
   #plt.close()

   capebins = np.arange(0, 501, 25)
   #umf_b_cape = bin_umf_by_var1(umf500, ds['CAPE'], capebins, ds['time'].size)
   #plt.plot(capebins, umf_b_cape)
   #plt.savefig('umf_binned_cape.png')
   #plt.close()

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

#not suited to work when bins don't cover the full range
def bin_umf_by_var1(umf, var1, bins, ntime):
   binidx = np.digitize(var1.data.ravel(), bins)
   return np.bincount(binidx, weights=umf.data.ravel(), minlength=bins.size) / ntime

def bin_umf_2d(umf, varx, vary, binsx, binsy, ntime):
   umfsum, xedges, yedges, binnum = binned_statistic_2d(x=varx.data.ravel(), y=vary.data.ravel(),\
                 values=umf.data.ravel(), statistic='sum', bins=[binsx, binsy])
   return umfsum / ntime, xedges, yedges, binnum

if __name__ == '__main__':
   main()
