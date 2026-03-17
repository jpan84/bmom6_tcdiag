import xarray as xr
import numpy as np
import os
import sys
from scipy.stats import binned_statistic
from datetime import timedelta as tdel
import matplotlib.pyplot as plt

ZMSST = 'zm_sst.nc' #full time series of zonal-mean SST from zm_sst_norm.py (case, time, yh)
pths = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm/hist', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist', '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist']
#pths = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm/hist', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/atm/hist']
dtstr = '*.h1i.0005-0[5-8]-*-*.nc'

g = 9.81
a = 6.371e6
umf500_f = lambda ds: (ds['OMEGA500'] < 0) * (-ds['OMEGA500'] * ds['area'] * a**2 / g)

def main():
   zmds = xr.open_dataset(ZMSST)
   maxlat = zmds['tos'].idxmax('yh')
   print(maxlat)
   sgnlat = np.sign(maxlat)
   prvsgn = sgnlat.shift(time=1)
   sgnflp = prvsgn == -sgnlat
   flip_dates = [sgnflp.sel(case=cs)['time'][sgnflp.sel(case=cs) == True] for cs in zmds['case']]

   print(flip_dates)
   #TODO: pick 1st occurrence within a ~7-day window

   test_dates = [da[0] for da in flip_dates]
   print(test_dates)

   dss = [compute_rel_SST(xr.open_mfdataset(os.path.join(pt, dtstr))) for pt in pths]
   dss = [ds.assign_coords(trel=('time', (ds['time'] - test_dates[ii]).data)) for ii, ds in enumerate(dss)]
   print(dss)

   histos = [compute_histo_tser(ds, thevarf=umf500_f)[0].assign_coords(trel=('time', (ds['time'] - test_dates[ii]).data)).expand_dims(case=[zmds['case'].isel(case=ii).data]) for ii, ds in enumerate(dss)]
   outda = xr.concat(histos, dim='case')
   outds = xr.Dataset(data_vars=dict(areasr=outda))
   outds.to_netcdf('0302_test_mf_vars/UMF500_SSTr_testhov.nc')

   print('hov complete')


def main_plot():
   testhov = '0302_test_mf_vars/areasr_SSTr_testhov.nc'
   ds = xr.open_dataset(testhov)
   #print(ds['trel'].isel(case=0).values)
   #print(ds.sel(trel=slice(tdel(days=-15), tdel(days=60))))

   ns2d = 1 / 1e9 / 86400
   selt = [ds.sel(case=cs) for cs in ds['case']]
   selt = [st.where((st['trel'].astype(float) * ns2d >= -15) & (st['trel'].astype(float) * ns2d <= 60), drop=True) for st in selt]
   selt = [st.swap_dims(time='trel') for st in selt]
   #print(selt[0])
   selt = [st - selt[1] if ii != 1 else st for ii, st in enumerate(selt)]  
   print(selt[0])

   plt.rcParams['figure.figsize'] = (12, 5)
   fig, axes = plt.subplots(1, 3, sharex=True, sharey=True)
   for ii, ax in enumerate(axes):
      csf = ax.contourf(ds['SSTr'], selt[ii]['trel'], selt[ii]['areasr'], cmap='bwr', levels=np.arange(-.3, .31, .05))
      plt.colorbar(csf, ax=ax)
   plt.show()

def compute_rel_SST(ds, troplat=30.):
   ds_trop = ds.isel(ncol=((ds['lat'] > -troplat) & (ds['lat'] < troplat)).all(dim='time').compute())
   meansst = (ds_trop['SST'] * ds_trop['area']).sum(dim='ncol') / ds_trop['area'].sum(dim='ncol')
   #print(meansst.compute())
   return ds.assign(SSTr=ds['SST'] - meansst)

def compute_histo_tser(ds, innerlat=5., outerlat=35., hemi='warm', thevarf=lambda ds: ds['area'], xvarf=lambda ds: ds['SSTr'],\
                     xbins=np.arange(-7, 5, 0.25), xnm='SSTr'):
   #print(ds.data_vars)
   selds = ds.isel(ncol=((ds['lat'] > innerlat) & (ds['lat'] < outerlat)).all(dim='time').compute()) #TODO: adapt to allow either hemi

   hist_snaps = []
   x_e = None
   for tt in selds['time']:
      print('Processing time', str(tt.data))
      binned, x_edges, _ = bin_stat(thevarf(selds).sel(time=tt), xvarf(selds).sel(time=tt),\
                                          xbins) #Divide time size by 2 because I've filtered out half of the year
      #print(umf_b_2d.shape, mse_edges.shape, sst_edges.shape)
      if x_e is None:
         x_e = x_edges

      hist_snaps.append(np.nan_to_num(binned, nan=0.0))

   #print(len(hist_snaps), hist_snaps[0].shape)
   hist_snaps = np.stack(hist_snaps, axis=0)
   x_c = (x_e[1:] + x_e[:-1]) / 2
   x_d = x_e[1:] - x_e[:-1]

   return xr.DataArray(hist_snaps, dims=['time', xnm], coords=[ds['time'], x_c]), x_d

def sel_unstruct_tropics(ds, innerlat=5., outerlat=30., hemi='warm'):
   szn = ds['time'].dt.season
   if hemi == 'warm':
      nh = ds.sel(time=szn.isin(['JJA', 'SON'])).isel(ncol=((ds['lat'] > innerlat) & (ds['lat'] < outerlat)).all(dim='time').compute())
      sh = ds.sel(time=szn.isin(['DJF', 'MAM'])).isel(ncol=((ds['lat'] < -innerlat) & (ds['lat'] > -outerlat)).all(dim='time').compute())
      return nh, sh
   if hemi == 'cool':
      nh = ds.sel(time=szn.isin(['DJF', 'MAM'])).isel(ncol=((ds['lat'] > innerlat) & (ds['lat'] < outerlat)).all(dim='time').compute())
      sh = ds.sel(time=szn.isin(['JJA', 'SON'])).isel(ncol=((ds['lat'] < -innerlat) & (ds['lat'] > -outerlat)).all(dim='time').compute())
      return nh, sh


def bin_stat(thevar, varx, binsx):
   summed, xedges, binnum = binned_statistic(x=varx.data.ravel(),\
                 values=thevar.data.ravel(), statistic='sum', bins=binsx)
   return summed, xedges, binnum

if __name__ == '__main__':
   if len(sys.argv) > 1 and sys.argv[1] == 'plot':
      main_plot()
   else:
      main()
