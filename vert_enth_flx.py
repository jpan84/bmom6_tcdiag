import xarray as xr
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from sznl_funcs import monthly2sznl, stack_hemi_sznl

DIRS = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm/hist_regrid_0.25x0.25_onpres/',\
        '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_regrid_0.25x0.25_onpres/',\
        '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist_regrid_0.25x0.25_onpres/']
ALIS = ['UNSEED', 'CTRL', 'MSEED']
DIRO = './w_enth_flx/'

FILI = 'ymonmean.nc' #'cdo_ann_means.nc' #'*h0a*.nc'
FILO = os.path.join(DIRO, 'w_enth_dT.nc')

PNM = 'plev'
cp = 1005
Rd = 287

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   print('Opening datasets...')
   dss = [xr.open_mfdataset(os.path.join(dr, FILI)).expand_dims(case=[ALIS[ii]]) for ii, dr in enumerate(DIRS)]
   inds = xr.concat(dss, dim='case').chunk(dict(plev=-1)).groupby('time.month').mean(dim='time')
   #print(inds)
   #exit()

   print('Taking zonal means...')
   vzm = inds['V'].mean(dim='lon')
   wzm = inds['OMEGA'].mean(dim='lon')
   tzm = inds['T'].mean(dim='lon')
   vt = inds['VT'].mean(dim='lon')
   wt = inds['OMEGAT'].mean(dim='lon')

   print('Reynolds decomposing')
   wt_mmc = wzm * tzm
   wt_tre = wt - wt_mmc #no stationary eddies
   vt_tre = vt - vzm * tzm
   
   print('Computing T tendencies...')
   dT_wadv_mmc = -wzm * tzm.differentiate(PNM, edge_order=2)
   dT_wwrk_mmc = Rd * wzm * tzm / cp / inds[PNM]
   dT_wflx_tre = -wt_tre.differentiate(PNM, edge_order=2)
   dT_wwrk_tre = Rd * wt_tre / cp / inds[PNM]

   print('Saving to .nc...')
   dvs = dict(dT_wadv_mmc=dT_wadv_mmc, dT_wwrk_mmc=dT_wwrk_mmc, dT_wflx_tre=dT_wflx_tre, dT_wwrk_tre=dT_wwrk_tre, OMEGAT_tre=wt_tre, VT_tre=vt_tre)
   outds = xr.Dataset(data_vars=dvs)
   outds.to_netcdf(FILO)

   print(sys.argv[0], 'done computing.')

def main_plot():
   ds = xr.open_dataset(FILO)

   dT_mmc = ds['dT_wadv_mmc'] + ds['dT_wwrk_mmc']
   dT_tre = ds['dT_wflx_tre'] + ds['dT_wflx_tre']

   m2s = lambda da: monthly2sznl(da)
   shs = lambda da: stack_hemi_sznl(da, antisym=False, latnm='lat')

   dT_mmc_sznl = shs(m2s(dT_mmc))
   dT_tre_sznl = shs(m2s(dT_tre))

   #plt.contourf(ds['lat'], ds['plev'], dT_tre_sznl.sel(season='SON', case='CTRL'), cmap='bwr', levels=np.arange(-8e-5, 8.1e-5, 1e-5))
   #plt.contourf(ds['lat'], ds['plev'], dT_tre_sznl.sel(season='SON', case='UNSEED') - dT_tre_sznl.sel(season='SON', case='CTRL'), cmap='RdBu_r', levels=np.arange(-4e-6, 4.1e-6, 5e-7))
   #plt.contourf(ds['lat'], ds['plev'], shs(m2s(-ds['OMEGAT_tre'])).sel(season='SON', case='CTRL'), cmap='bwr', levels=np.arange(-0.25, 0.26, .05))
   plt.contourf(ds['lat'], ds['plev'], shs(m2s(-ds['OMEGAT_tre'])).sel(season='SON', case='UNSEED') - shs(m2s(-ds['OMEGAT_tre'])).sel(season='SON', case='CTRL'), cmap='RdBu_r', levels=np.arange(-0.01, 0.011, .002))
   plt.colorbar()
   plt.xlim(-30, 30)
   plt.ylim(1e5, 1e4)
   plt.show()

def main_vT():
   ds = xr.open_dataset(FILO)
   vt850_sznl = stack_hemi_sznl(monthly2sznl(ds['VT_tre'].interp(plev=8.5e4)), antisym=True, latnm='lat')
   #plt.plot(ds['lat'], vt850_sznl.sel(season='SON', case='CTRL'))
   plt.plot(ds['lat'], vt850_sznl.sel(season='SON', case='UNSEED') - vt850_sznl.sel(season='SON', case='CTRL'))
   plt.show()

if __name__ == '__main__':
   if sys.argv[1] == 'compute':
      main()
   if sys.argv[1] == 'plot':
      main_plot()
   if sys.argv[1] == 'vT':
      main_vT()
