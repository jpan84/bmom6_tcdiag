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
FILIS = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/atm/hist/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s.cam.h1i.0012-*-0*-*.nc'
ALIS = ['250415_unseed', '250417_ctrl', '250416_seed1x1']
CASES = ['UNSEED', 'CTRL', 'SEED']
#FILO = 'umf500_0012dd0x_histo_mse850_SST.pkl'

SSTbins = np.concatenate((np.arange(295., 306.), np.arange(306, 309, 0.25), np.arange(309, 313)))

#TODO: allow flexible latitude bins?
#TODO: allow diffing cases
#TODO: split seasons
#TODO: check why NH and SH selections differ in number of cols
def main():
   #print(sys.argv)
   #ds = xr.open_mfdataset(CASE1)


   mtout, stout = None, None
   m_d, s_d = None, None
   for ii, al in enumerate(ALIS):
      print('Working on', al)
      ds = xr.open_mfdataset(FILIS % (al, al))
      mtda, m_d, s_d = compute_umf_hist(ds.drop_vars(['time_bounds', 'time_written', 'date_written']), hemi='warm', sstbins=SSTbins, innerlat=7., outerlat=20.)
      stda, _, _ = compute_umf_hist(ds.drop_vars(['time_bounds', 'time_written', 'date_written']), hemi='warm', sstbins=SSTbins, innerlat=20., outerlat=33.)
      mtda, stda = mtda.expand_dims(case=[CASES[ii]]), stda.expand_dims(case=[CASES[ii]])

      if mtout is None:
         mtout, stout = mtda, stda
      else:
         mtout, stout = xr.concat([mtout, mtda], dim='case'), xr.concat([stout, stda], dim='case')
      #topkl.append(compute_umf_hist(ds, hemi='cool', capebins=SSTbins))

   mtds = xr.Dataset(data_vars=dict(UMF500=mtout, mwidth=m_d, swidth=s_d))
   stds = xr.Dataset(data_vars=dict(UMF500=stout, mwidth=m_d, swidth=s_d))

   mtds.to_netcdf('umf500_0012dd0x_07-20warm.nc')
   stds.to_netcdf('umf500_0012dd0x_20-33warm.nc')

   print(sys.argv[0], 'done')



def compute_umf_hist(ds, innerlat=5., outerlat=30., hemi='warm', msebins=np.arange(2.8e5, 4.01e5, 5e3), sstbins=np.arange(295, 315)): #try SST, TMQ
   #print(ds.data_vars)
   nhsel, shsel = sel_unstruct_tropics(ds, innerlat=innerlat, outerlat=outerlat, hemi=hemi)
   selds = xr.concat((nhsel, shsel), 'ncol')
   #print(selds)
   #print(selds.time)

   mse850 = g * selds['Z850'] + cp * selds['T850'] + lv * selds['Q850']
   area = selds['area'] * a**2
   umf500 = (selds['OMEGA500'] < 0) * (-selds['OMEGA500'] * area / g)

   hist_snaps = []
   m_e, s_e = None, None
   for tt in selds['time']:
      print('Processing time', str(tt.data))
      umf_b_2d, mse_edges, sst_edges, _ = bin_umf_2d(umf500.sel(time=tt), mse850.sel(time=tt), selds['SST'].sel(time=tt),\
                                          msebins, sstbins, selds['time'].size / 2) #Divide time size by 2 because I've filtered out half of the year
      #print(umf_b_2d.shape, mse_edges.shape, sst_edges.shape)
      if m_e is None and s_e is None:
         m_e, s_e = mse_edges, sst_edges

      hist_snaps.append(np.nan_to_num(umf_b_2d, nan=0.0))

   #print(len(hist_snaps), hist_snaps[0].shape)
   hist_snaps = np.stack(hist_snaps, axis=0)
   m_c, s_c = (m_e[1:] + m_e[:-1]) / 2, (s_e[1:] + s_e[:-1]) / 2
   m_d, s_d = m_e[1:] - m_e[:-1], s_e[1:] - s_e[:-1]

   return xr.DataArray(hist_snaps, dims=['time', 'MSE850', 'SST'], coords=[selds['time'], m_c, s_c]), m_d, s_d

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
