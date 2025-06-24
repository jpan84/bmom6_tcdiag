#python3 bin_mass_flux.py 0 cool
#python3 bin_mass_flux.py 1 warm
import sys
import numpy as np
from scipy.stats import binned_statistic_2d
#import uxarray as ux
import xarray as xr
import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors

g = 9.81
cp = 1004
lv = 2500840
a = 6.371e6

#camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

#DO_DIFF = bool(int(sys.argv[1])) #True
#hem = sys.argv[2]
#alias2 = sys.argv[3] if DO_DIFF else None

#CASE1 = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_0010_h1i/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.*.nc'
#CASE2 = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/atm/hist_0010_h1i/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s.cam.h1i.*.nc'\
#         % (alias2, alias2)
FILIS = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/atm/hist_0010_h1i/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s.cam.h1i.*.nc'
ALIS = ['250415_unseed', '250417_ctrl', '250416_seed1x1']
FILO = 'umf500_0010_histo_mse850_cape.pkl'

#TODO: allow flexible latitude bins?
#TODO: allow diffing cases
#TODO: split seasons
#TODO: check why NH and SH selections differ in number of cols
def main():
   #print(sys.argv)
   #ds = xr.open_mfdataset(CASE1)

   topkl = [tuple(ALIS), ('warm', 'cool')]

   for al in ALIS:
      print('Working on', al)
      ds = xr.open_mfdataset(FILIS % (al, al))
      topkl.append(compute_umf_hist(ds, hemi='warm'))
      topkl.append(compute_umf_hist(ds, hemi='cool'))

   print('Pickling')
   with open(FILO, 'wb') as fl:
      pickle.dump(topkl, fl)

   print(sys.argv[0], 'done')

   exit()
   ##################################
   print('Computing histogram...')
   umf_b_2d, mse_edges, cape_edges = compute_umf_hist(ds, hemi=hem)
   if DO_DIFF:
      ds2 = xr.open_mfdataset(CASE2)
      umf_b_2d_CASE2, _, _ = compute_umf_hist(ds2, hemi=hem)
      umf_b_2d = umf_b_2d_CASE2 - umf_b_2d

   print('Plotting...')
   plt.rc('font', size=16)
   #print(mse_edges, cape_edges, umf_b_2d.shape)
   vminmax = dict(vmin=0, vmax=1e11)
   if DO_DIFF:
      vminmax = dict()
   plt.pcolormesh(mse_edges, cape_edges, umf_b_2d.T, shading='flat', cmap='bwr' if DO_DIFF else 'viridis', **vminmax)#, norm=colors.LogNorm())
   plt.xlim(3.0e5, 3.8e5)
   plt.ylim(0, 400)
   plt.xlabel('MSE 850 [J kg$^{-1}$]')
   plt.ylabel('CAPE [J kg$^{-1}$]')
   plt.title(f'500 hPa UMF [kg s$^{{-1}}$]\n{umf_b_2d.sum():.3e}')
   plt.colorbar()
   plt.savefig('umf_binned_2d_%s%s.png' % (hem, '_%s_diff' % alias2 if DO_DIFF else ''), bbox_inches='tight')
   plt.close()

   print('UMFBIN done.')

def compute_umf_hist(ds, outerlat=30, hemi='warm', msebins=np.arange(2.8e5, 4.01e5, 5e3), capebins=np.arange(0, 1001, 25)):
   nhsel, shsel = sel_unstruct_tropics(ds, outerlat=outerlat, hemi=hemi)
   selds = xr.concat((nhsel, shsel), 'ncol')

   mse850 = g * selds['Z850'] + cp * selds['T850'] + lv * selds['Q850']
   area = selds['area'] * a**2
   umf500 = (selds['OMEGA500'] < 0) * (-selds['OMEGA500'] * area / g)

   umf_b_2d, mse_edges, cape_edges, _ = bin_umf_2d(umf500, mse850, selds['CAPE'], msebins, capebins, selds['time'].size / 2) #Divide time size by 2 because I've filtered out half of the year
   return np.nan_to_num(umf_b_2d, nan=0.0), mse_edges, cape_edges


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
