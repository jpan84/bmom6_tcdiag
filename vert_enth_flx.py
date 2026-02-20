import xarray as xr
import os
import sys
import matplotlib.pyplot as plt

DIRS = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_regrid_0.25x0.25_onpres/',\
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
   print('Opening datasets...')
   dss = [xr.open_mfdataset(os.path.join(dr, FILI)).expand_dims(case=[ALIS[ii]]) for ii, dr in enumerate(DIRS)]
   inds = xr.concat(dss, dim='case')

   print('Taking zonal means...')
   wzm = inds['OMEGA'].mean(dim='lon')
   tzm = inds['T'].mean(dim='lon')
   wt = inds['OMEGAT'].mean(dim='lon')

   print('Reynolds decomposing')
   wt_mmc = wzm * tzm
   wt_tre = wt - wt_mmc #no stationary eddies
   
   print('Computing T tendencies...')
   dT_wadv_mmc = -wzm * tzm.differentiate(PNM, edge_order=2)
   dT_wwrk_mmc = Rd * wzm * tzm / cp / inds[PNM]
   dT_wflx_tre = -wt_tre.differentiate(PNM, edge_order=2)
   dT_wwrk_tre = Rd * wt_tre / cp / inds[PNM]

   print('Saving to .nc...')
   dvs = dict(dT_wadv_mmc=dT_wadv_mmc, dT_wwrk_mmc=dT_wwrk_mmc, dT_wflx_tre=dT_wflx_tre, dT_wwrk_tre=dT_wwrk_tre)
   outds = xr.Dataset(data_vars=dvs)
   outds.to_netcdf(FILO)

   print(sys.argv[0], 'done computing.')

def main_plot():
   return

if __name__ == '__main__':
   if sys.argv[1] == 'compute':
      main()
   if sys.argv[1] == 'plot':
      main_plot()
