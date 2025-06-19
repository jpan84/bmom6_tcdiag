import numpy as np
import uxarray as ux
import matplotlib.pyplot as plt

g = 9.81
cp = 1004
lv = 2500840
a = 6.371e6

testfile = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_0010_h1i/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.0010-07-0*.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

def main():
   ds = ux.open_mfdataset(camgrid, testfile)
   #print(ds['area'].min().values, ds['area'].max().values) 

   mse850 = g * ds['Z850'] + cp * ds['T850'] + lv * ds['Q850']
   plt.hist(mse850.data)
   plt.savefig('mse850_hist.png')
   plt.close()

   area = ds['area'] * a**2
   umf500 = (ds['OMEGA500'] < 0) * (-ds['OMEGA500'] * area / g)

   msebins = np.arange(2.5e5, 4e5, 2.5e4)
   ###msebinidx = np.digitize(np.ravel(mse850.data), msebins)

   ####serial way
   ######umf_binned = np.empty_like(msebins)
   ######for ii in range(umf_binned.size):
   ######   inbin = np.where(msebinidx == ii)[0]
   ######   umf_binned[ii] = umf500.data.ravel()[inbin].sum()

   ####vectorized
   ###umf_binned = np.bincount(msebinidx, weights=umf500.data.ravel(), minlength=msebins.size) / ds['time'].size

   umf_b_h850 = bin_umf_by_var1(umf500, mse850, msebins, ds['time'].size)
   plt.plot(msebins, umf_b_h850)
   plt.savefig('umf_binned_mse850.png')
   plt.close()

   capebins = np.arange(-400, 1001, 200)
   umf_b_cape = bin_umf_by_var1(umf500, ds['CAPE'], capebins, ds['time'].size)
   plt.plot(capebins, umf_b_cape)
   plt.savefig('umf_binned_cape.png')
   plt.close()


def bin_umf_by_var1(umf, var1, bins, ntime):
   binidx = np.digitize(var1.data.ravel(), bins)
   return np.bincount(binidx, weights=umf.data.ravel(), minlength=bins.size) / ntime

if __name__ == '__main__':
   main()
