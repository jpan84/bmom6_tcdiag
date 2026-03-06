#python3 bin_mass_flux.py 0 cool
#python3 bin_mass_flux.py 1 warm
import sys
import os
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
FILIS = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/atm/hist/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s.cam.h1i.0007-*-1*-*.nc'
ALIS = ['250415_unseed', '250417_ctrl', '251229_seedmatch']#'250416_seed1x1']
pths = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm/hist', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist', '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist']
dtstr = '*.h1i.0009-*-01-*.nc' #'*.h1i.000[6-7]-*.nc'
CASES = ['UNSEED', 'CTRL', 'MSEED']
VARO = 'UMF500'
ILAT, OLAT = 5, 35

mse850_f = lambda ds: g * ds['Z850'] + cp * ds['T850'] + lv * ds['Q850']
spdbot_f = lambda ds: np.sqrt(ds['UBOT']**2 + ds['VBOT']**2)
sst_f = lambda ds: ds['SST']
lat_f = lambda ds: ds['lat']
umf500_f = lambda ds: (ds['OMEGA500'] < 0) * (-ds['OMEGA500'] * ds['area'] * a**2 / g)
umf850_f = lambda ds: (ds['OMEGA850'] < 0) * (-ds['OMEGA850'] * ds['area'] * a**2 / g) #* (ds['PRECT'] > 1e-8)
dmf500_f = lambda ds: (ds['OMEGA500'] > 0) * (-ds['OMEGA500'] * ds['area'] * a**2 / g)

SSTbins = np.concatenate((np.arange(295., 306.), np.arange(306, 309, 0.25), np.arange(309, 313)))
OM850b = np.concatenate((np.arange(-10, -2, 1), np.arange(-2, -.5, .25), np.arange(-.5, -.1, .05), np.arange(-.1, .01, .01)))
#OM850b = np.arange(-40, 0.25, 1)
OM500b = np.concatenate((np.arange(-20, -4, 2), np.arange(-4, -.5, .25), np.arange(-.5, -.2, .1), np.arange(-.2, 0, .02), np.arange(0, .51, .1)))
#OM500b = np.arange(-80, 10, 2)
FLUTb = np.concatenate((np.arange(80, 180, 5), np.arange(180, 361, 5)))
latb_500u = np.concatenate((np.arange(5, 8, 0.5), np.arange(8, 14, 0.25), np.arange(14, 20, 0.5), np.arange(20, 36, 1)))
latb_500d = np.arange(5, 36, 1)
UBOTb = np.arange(0, 41, 2)

hist_kw = dict(hemi='warm', xnm='FLUT', ynm='OM500', innerlat=ILAT, outerlat=OLAT, thevarf=umf850_f, xvarf=lambda x: x['FLUT'], yvarf=lambda x: x['OMEGA500'], xbins=FLUTb, ybins=OM500b)
hist_kw = dict(hemi='warm', xnm='OM850', ynm='OM500', innerlat=ILAT, outerlat=OLAT, thevarf=umf850_f, xvarf=lambda x: x['OMEGA850'], yvarf=lambda x: x['OMEGA500'], xbins=OM850b, ybins=OM500b)
hist_kw = dict(hemi='warm', xnm='lat', ynm='SST', innerlat=ILAT, outerlat=OLAT, thevarf=umf500_f, xvarf=lat_f, yvarf=sst_f, xbins=latb_500u, ybins=SSTbins)
hist_kw = dict(hemi='warm', xnm='lat', ynm='SST', innerlat=ILAT, outerlat=OLAT, thevarf=dmf500_f, xvarf=lat_f, yvarf=sst_f, xbins=latb_500d, ybins=SSTbins)
hist_kw = dict(hemi='warm', xnm='lat', ynm='UBOT', innerlat=ILAT, outerlat=OLAT, thevarf=umf500_f, xvarf=lat_f, yvarf=spdbot_f, xbins=latb_500u, ybins=UBOTb)

#TODO: check why NH and SH selections differ in number of cols
def main():
   global SSTbins

   outda = None
   x_d, y_d = None, None
   for ii, al in enumerate(ALIS):
      print('Working on', al)
      #ds = xr.open_mfdataset(FILIS % (al, al))
      ds = xr.open_mfdataset(os.path.join(pths[ii], dtstr))
      #da, x_d, y_d = compute_umf_hist(ds.drop_vars(['time_bounds', 'time_written', 'date_written']), hemi='warm', ybins=SSTbins, innerlat=7., outerlat=20.)
      #da, x_d, y_d = compute_umf_hist(ds.drop_vars(['time_bounds', 'time_written', 'date_written']), hemi='warm', innerlat=ILAT, outerlat=OLAT,\
      #                                  thevarf=umf850_f, xvarf=lambda x: x['OMEGA850'], yvarf=lambda x: x['OMEGA500'], xbins=OM850b, ybins=OM500b, xnm='OM850', ynm='OM500')
      da, x_d, y_d = compute_umf_hist(ds.drop_vars(['time_bounds', 'time_written', 'date_written']), **hist_kw)
      da = da.expand_dims(case=[CASES[ii]])

      if outda is None:
         outda = da
      else:
         outda = xr.concat([outda, da], dim='case')

   outds = xr.Dataset(data_vars={VARO: outda, 'xwidth': x_d, 'ywidth': y_d})
   outds.to_netcdf('0302_test_mf_vars/%s_%s_%s_yy09dd01_%02d-%02d%s_allprec.nc' % (VARO, hist_kw['xnm'], hist_kw['ynm'], ILAT, OLAT, hist_kw['hemi']))

   print(sys.argv[0], 'done')



def compute_umf_hist(ds, innerlat=5., outerlat=30., hemi='warm', thevarf=umf500_f, xvarf=mse850_f, yvarf=sst_f,\
                     xbins=np.arange(2.8e5, 4.01e5, 5e3), ybins=np.arange(295, 315), xnm='MSE850', ynm='SST'):
   #print(ds.data_vars)
   nhsel, shsel = sel_unstruct_tropics(ds, innerlat=innerlat, outerlat=outerlat, hemi=hemi)
   selds = xr.concat((nhsel, shsel), 'ncol')

   hist_snaps = []
   m_e, s_e = None, None
   for tt in selds['time']:
      print('Processing time', str(tt.data))
      umf_b_2d, mse_edges, sst_edges, _ = bin_umf_2d(thevarf(selds).sel(time=tt), xvarf(selds).sel(time=tt), yvarf(selds).sel(time=tt),\
                                          xbins, ybins, selds['time'].size / 2) #Divide time size by 2 because I've filtered out half of the year
      #print(umf_b_2d.shape, mse_edges.shape, sst_edges.shape)
      if m_e is None and s_e is None:
         m_e, s_e = mse_edges, sst_edges

      hist_snaps.append(np.nan_to_num(umf_b_2d, nan=0.0))

   #print(len(hist_snaps), hist_snaps[0].shape)
   hist_snaps = np.stack(hist_snaps, axis=0)
   m_c, s_c = (m_e[1:] + m_e[:-1]) / 2, (s_e[1:] + s_e[:-1]) / 2
   m_d, s_d = m_e[1:] - m_e[:-1], s_e[1:] - s_e[:-1]

   return xr.DataArray(hist_snaps, dims=['time', xnm, ynm], coords=[selds['time'], m_c, s_c]), m_d, s_d

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


def bin_umf_2d(umf, varx, vary, binsx, binsy, ntime):
   umfsum, xedges, yedges, binnum = binned_statistic_2d(x=varx.data.ravel(), y=vary.data.ravel(),\
                 values=umf.data.ravel(), statistic='sum', bins=[binsx, binsy])
   return umfsum / ntime, xedges, yedges, binnum

if __name__ == '__main__':
   main()
