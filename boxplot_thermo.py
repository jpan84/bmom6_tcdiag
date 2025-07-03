#python3 bin_mass_flux.py 0 cool
#python3 bin_mass_flux.py 1 warm
import os
import sys
import numpy as np
import dask
import operator
#import uxarray as ux
import xarray as xr
import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import metpy
from metpy.plots import SkewT

g = 9.81
cp = 1004
lv = 2500840
a = 6.371e6

#camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
DZTHRESH = 10 #mm/d threshold to separate non/precipitating regimes
OPs = {'<': operator.lt, '<=': operator.le, '>': operator.gt, '>=': operator.ge}
VARPLEVS = dict(OMEGA=[500, 850], T=[200, 500, 850, 1010], Z=[300, 500, 850], U=[200, 500, 850])
OFF = {200: 10, 500: 25, 850: 40, 1010: 50}

HEM = 'cool'
MYOP = '>='
FILIS = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/atm/hist_0010_h1i/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s.cam.h1i.0010-*-01-*.nc'
ALIS = ['250415_unseed', '250417_ctrl', '250416_seed1x1']
PLCLRS = ['blue', 'green', 'red']
DIRO = './thermo_boxplots'

#TODO: check why NH and SH selections differ in number of cols
def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)
   dss = [xr.open_mfdataset(FILIS % (al, al)).rename(name_dict=dict(TBOT='T1010')) for al in ALIS]

   topkl = [tuple(ALIS), ('warm', 'cool')]

   for var in VARPLEVS:
      skew = None
      if var == 'T':
         skew = metpy.plots.SkewT()
         skew.plot_moist_adiabats()
      for jj, pl in enumerate(VARPLEVS[var]):
         for ii, al in enumerate(ALIS):
            print('Working on', al, var, pl)
            if not var == 'T': continue #TODO
            ds = dss[ii]
            quartiles = get_weighted_quartiles_filtered(ds, 'PRECT', MYOP, DZTHRESH / 8.64e4 / 1e3, '%s%d' % (var, pl), hemi=HEM)
            #plt.boxplot(quartiles)
            if var == 'T':
               offset = OFF[pl] #25 if pl >= 500 else 10
               plt.gca().scatter(quartiles - 273.15, np.zeros_like(quartiles) + pl + (ii - 1) * offset, c=PLCLRS[ii], marker='|')
               #quartiles = quartiles * (1000. / pl)**(287 / 1004)
            else:
               plt.scatter(quartiles, np.zeros_like(quartiles) + pl + (ii - 1) * 15, c=PLCLRS[ii], marker='|')

      plt.xlabel(var)
      #if var == 'T':
      #   plt.xlabel('Pot temp [K]')
      plt.ylabel('Pressure [hPa]')
      plt.yticks(VARPLEVS[var])
      plt.gca().invert_yaxis()
      if var == 'T':
         plt.gca().set_ylim(1100, 180)
         plt.gca().set_xlim(-20, 40)
      plt.savefig(os.path.join(DIRO, '%s_%s%d_%s.png' % (var, OPs[MYOP].__name__, DZTHRESH, HEM)))
      plt.close()

      continue
      topkl.append(compute_umf_hist(ds, hemi='warm'))
      topkl.append(compute_umf_hist(ds, hemi='cool'))

   exit()
   print('Pickling')
   with open(FILO, 'wb') as fl:
      pickle.dump(topkl, fl)

   print(sys.argv[0], 'done')

   exit()


#use one variable and threshold to filter datapoints
#get the area-weighted CDF/quartiles of any var in the filtered dataset
def get_weighted_quartiles_filtered(ds, varfilt, oper, thresh, cdfvar, wgts=None, outerlat=30., hemi='warm'):
   nhsel, shsel = sel_unstruct_tropics(ds, outerlat=outerlat, hemi=hemi)
   selds = xr.concat((nhsel, shsel), 'ncol')
   #print(selds)

   mask = OPs[oper](selds[varfilt], thresh).data.reshape(-1).compute()
   wgts = selds['area'].data.reshape(-1)[mask].compute()
   var1d = selds[cdfvar].data.reshape(-1)[mask].compute()

   #mask = OPs[oper](selds[varfilt], thresh).compute()
   #wgts = dask.array.from_array(selds['area'].where(mask, drop=False))
   #var1d = dask.array.from_array(selds[cdfvar].where(mask, drop=False))

   return np.nanquantile(var1d, np.arange(.05, .96, 0.1), weights=wgts, method='inverted_cdf') #compute_wgted_quartiles(var1d, wgts)

###def compute_wgted_quartiles(var, wgts):
###   sorter = np.argsort(var)
###   varsort = var[sorter]
###   wgtssort = wgts[sorter]
###   quantiles = np.cumsum(wgts) / sum(wgts)
###   return np.interp(np.arange(0, 1.01, 0.25), quantiles, varsort)

#TODO: check why the 2 hemis return different numbers of points
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


if __name__ == '__main__':
   main()
