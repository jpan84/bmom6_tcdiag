import xarray as xr
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sznl_funcs import monthly2sznl, stack_hemi_sznl

DIRS = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm/hist_regrid_0.25x0.25_onpres/',\
        '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_regrid_0.25x0.25_onpres/',\
        '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist_regrid_0.25x0.25_onpres/']
ALIS = ['UNSEED', 'CTRL', 'MSEED']
DIRO = './w_enth_flx/'

FILI = 'ymonmean.nc' #'cdo_ann_means.nc' #'*h0a*.nc'
FILO = os.path.join(DIRO, 'w_enth_dT_vflxs_diab.nc')

PNM = 'plev'
cp = 1005
Rd = 287

YSCL = lambda lat: np.sin(np.deg2rad(lat))

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

   print('Including v fluxes...')
   vu = inds['VU'].mean(dim='lon')
   vq = inds['VQ'].mean(dim='lon')
   uzm = inds['U'].mean(dim='lon')
   qzm = inds['Q'].mean(dim='lon')

   print('Including diab...')
   dT_LH = inds['DTCOND'].mean(dim='lon')
   dT_LW = inds['QRL'].mean(dim='lon')
   dT_SW = inds['QRS'].mean(dim='lon')

   print('Reynolds decomposing...')
   wt_mmc = wzm * tzm
   wt_tre = wt - wt_mmc #no stationary eddies
   vt_tre = vt - vzm * tzm
   vu_tre = vu - vzm * uzm
   vq_tre = vq - vzm * qzm
   
   print('Computing T tendencies...')
   dT_wadv_mmc = -wzm * tzm.differentiate(PNM, edge_order=2)
   dT_wwrk_mmc = Rd * wzm * tzm / cp / inds[PNM]
   dT_wflx_tre = -wt_tre.differentiate(PNM, edge_order=2)
   dT_wwrk_tre = Rd * wt_tre / cp / inds[PNM]

   print('Saving to .nc...')
   dvs = dict(dT_wadv_mmc=dT_wadv_mmc, dT_wwrk_mmc=dT_wwrk_mmc, dT_wflx_tre=dT_wflx_tre, dT_wwrk_tre=dT_wwrk_tre, OMEGAT_tre=wt_tre,\
              VT_tre=vt_tre, VU_tre=vu_tre, VQ_tre=vq_tre, dT_LH=dT_LH, dT_LW=dT_LW, dT_SW=dT_SW)
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
   #plt.contourf(ds['lat'], ds['plev'], shs(m2s(-ds['OMEGAT_tre'])).sel(season='SON', case='UNSEED') - shs(m2s(-ds['OMEGAT_tre'])).sel(season='SON', case='CTRL'), cmap='RdBu_r', levels=np.arange(-0.01, 0.011, .002))
   #plt.contourf(ds['lat'], ds['plev'], shs(m2s(ds['dT_LW'])).sel(season='SON', case='MSEED') - shs(m2s(ds['dT_LW'])).sel(season='SON', case='CTRL'), cmap='RdBu_r', levels=np.arange(-1e-6, 1.1e-6, 2e-7))
   #plt.contourf(ds['lat'], ds['plev'], shs(m2s(ds['dT_LH'])).sel(season='SON', case='UNSEED') - shs(m2s(ds['dT_LH'])).sel(season='SON', case='CTRL'), cmap='RdBu_r', levels=np.arange(-5e-6, 5.1e-6, 1e-6))

   plt.contourf(ds['lat'], ds['plev'], shs(m2s(-ds['OMEGAT_tre'])).sel(season='SON', case='MSEED') - shs(m2s(-ds['OMEGAT_tre'])).sel(season='SON', case='CTRL'), cmap='RdBu_r', levels=np.arange(-0.01, 0.011, .002))
   plt.colorbar()
   plt.contour(ds['lat'], ds['plev'], shs(m2s(ds['dT_LH'])).sel(season='SON', case='MSEED') - shs(m2s(ds['dT_LH'])).sel(season='SON', case='CTRL'), colors='black', levels=np.arange(-5e-6, 5.1e-6, 1e-6))
   plt.xlim(-30, 30)
   plt.ylim(1e5, 1e4)
   plt.yscale('log')
   plt.show()

def main_vT():
   ds = xr.open_dataset(FILO)
   usds = xr.open_dataset('./tcfields4mps_250415_unseed/means_tcs.nc')
   ctds = xr.open_dataset('./tcfields4mps_250417_ctrl/means_tcs.nc')
   vt850_net = lambda ds: ds['VT850_neg'] + ds['VT850_pos']
   pltszn = 'JJA'

   vt850_sznl = stack_hemi_sznl(monthly2sznl(ds['VT_tre'].interp(plev=8.5e4)), antisym=True, latnm='lat')
   toplt_tc = (vt850_net(usds) - vt850_net(ctds)).sel(season=pltszn)
   #plt.plot(ds['lat'], vt850_sznl.sel(season='SON', case='CTRL'))
   plt.plot(ds['lat'], vt850_sznl.sel(season=pltszn, case='UNSEED') - vt850_sznl.sel(season=pltszn, case='CTRL'))
   plt.plot(ctds['latitudes'], toplt_tc)
   plt.show()

#plot v eddy transports
def main_vX():
   ds = xr.open_dataset(FILO)
   ds = ds.assign(plev=ds['plev'] / 100)
   pltvars = ['VT_tre', 'VU_tre', 'VQ_tre']
   clevs_CTRL = dict(VU_tre=np.arange(-100, 101, 10), VT_tre=np.arange(-10, 11, 2), VQ_tre=np.arange(-1e-2, 1.1e-2, 1e-3))
   clevs_diff = dict(VU_tre=np.arange(-2.5, 2.6, 0.5), VT_tre=np.arange(-0.5, 0.51, 0.1), VQ_tre=np.arange(-5e-4, 5.1e-4, 1e-4))

   m2s = lambda da: monthly2sznl(da)
   shs = lambda da: stack_hemi_sznl(da, antisym=True, latnm='lat')  

   flds_sznl = [shs(m2s(ds[pv])) for pv in pltvars]

   plt.rc('font', size=20)
   plt.rcParams['figure.figsize'] = (30, 12)
   latlim = YSCL(50)
   subplot_kw = dict(xlim=(-latlim, latlim), ylim=(100, 1000), yscale='log')

   for vrk, fld in enumerate(flds_sznl):
      fld.loc[dict(case=['UNSEED', 'MSEED'])] -= fld.sel(case='CTRL')
      fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, subplot_kw=subplot_kw) #layout='constrained'
      axes[0][0].invert_yaxis()
   
      for csi, ali in enumerate(ALIS):
         for szj, szn in enumerate(fld['season']):
            ax, do_diff = axes[szj][csi], ali != 'CTRL'
            csf = ax.contourf(YSCL(fld['lat']), fld['plev'], fld.sel(case=ali, season=szn), norm=mcolors.CenteredNorm(), levels=clevs_diff[pltvars[vrk]] if do_diff else clevs_CTRL[pltvars[vrk]], cmap='bwr')
            plt.colorbar(csf)
            ax.set_xticks(YSCL(np.arange(-50, 51, 10)), np.arange(-50, 51, 10))
      plt.savefig(os.path.join(DIRO, pltvars[vrk]))
      plt.close()


if __name__ == '__main__':
   if sys.argv[1] == 'compute':
      main()
   if sys.argv[1] == 'plot':
      main_plot()
   if sys.argv[1] == 'vT':
      main_vT()
   if sys.argv[1] == 'vX':
      main_vX()
