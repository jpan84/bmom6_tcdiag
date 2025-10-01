#Step 1: obtain hybrid-level pressure thickness from hyb_delta_pressure.py
#Step 2: run cdo commands to obtain monthly norms of vertically integrated meridional fluxes. e.g.,
#   cdo selvar,VQ,VT,VU,Q,T,U,Z3,V -cat *h0a*.nc cat_merid.nc
#   setsid nohup cdo --debug ymonmean -vertsum -mul [ [ -selvar,VQ cat_merid.nc ] [ -selvar,dp3d dp_monthly.nc ] ] vqdp_monthly.nc &> vqdp_monthly.out
#   setsid nohup cdo --debug ymonmean -vertsum -mul [ [ -mul [ [ -selvar,V cat_merid.nc ] [ -selvar,Q cat_merid.nc ] ] ] [ -selvar,dp3d dp_monthly.nc ] ] vqdp_mmc_monthly.nc &> vqdp_mmc_monthly.out

import sys
import os
import math
import numpy as np
import xarray as xr
import uxarray as ux
from dask.diagnostics import ProgressBar
from sznl_funcs import monthly2sznl, stack_hemi_sznl
import matplotlib.pyplot as plt

MODE = sys.argv[1] #'compute', 'plot'

FILO = 'transports_3D_conserv.nc'
DIRO = 'merid_tprt_conserv/'

CAMGR = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
CAMDIR = '/glade/campaign/univ/upsu0032/jpan_aquaptc//%s/atm/hist/'
MOMH = '/glade/derecho/scratch/jpan/archive/%s/ocn/hist/*mom6.hm*[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'
CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
ALIASES = ['UNSEED', 'CTRL', 'SEED']

P0 = 1e5
g = 9.79764
OM = 7.2921e-5
a_e = 6.371e6
cp = 1005
lv = 2.501e6
LATLAB = np.array([-90., -60., -50., -40., -30., -20, -10, 0., 10., 20., 30., 40., 50., 60., 90.])

#inflds = [('Q', 'VQ')]#, ('T', 'VT'), ('U', 'VU'), ('Z3', None)]
tprtfiles = dict(q=('vqdp_monthly.nc', 'vqdp_mmc_monthly.nc'), T=('vTdp_monthly.nc', 'vTdp_mmc_monthly.nc'), u=('vudp_monthly.nc', 'vudp_mmc_monthly.nc'), Z=(None, 'vZdp_mmc_monthly.nc'))

def main():
   print('Processing total transports...')
   trtot = [{kk: ux.open_mfdataset(CAMGR, os.path.join(CAMDIR % cs, vv[1])) if vv[0] is None else ux.open_mfdataset(CAMGR, os.path.join(CAMDIR % cs, vv[0])) for kk, vv in tprtfiles.items()} for cs in CASES]
   #print([list(dd['q'].data_vars.values())[1] for dd in trtot])
   trtot = [{kk: list(vv.data_vars.values())[-1] for kk, vv in dd.items()} for dd in trtot] #pick out the data_var containing the transport field (at idx -1)
   #print(trtot)
   trtot = [ux.UxDataset(data_vars=dd).groupby('time.month').mean('time') for dd in trtot] #gather the different fields into a single Dataset for each case
   #print(trtot)
   #casecoords = xr.Coordinates(coords=dict(case=ALIASES))
   totds = xr.concat([ds.expand_dims(case=[ALIASES[ii]]) for ii, ds in enumerate(trtot)], dim='case')
   print(totds)

   '''
   trmmc = []

   # Loop 1: Iterate over each case (e.g., 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed')
   for cs in CASES:
      # This dictionary will hold the DataArrays for the current case (cs)
      case_data = {}
      
      # Loop 2: Iterate over each variable and its corresponding file paths
      # kk is the variable name (e.g., 'q')
      # vv is the tuple of file names (e.g., ('vqdp_monthly.nc', 'vqdp_mmc_monthly.nc'))
      for kk, vv in tprtfiles.items():
         # vv[1] is the second file name (e.g., 'vqdp_mmc_monthly.nc')
         file_name = vv[1]
         
         # Construct the full path
         full_path = os.path.join(CAMDIR % cs, file_name)
         
         print(f"Attempting to open: {full_path}") # Debugging line
         
         try:
            # Open the dataset for the current variable and case
            ds = ux.open_mfdataset(CAMGR, full_path)
            
            # Store the resulting Dataset in the case's dictionary
            case_data[kk] = ds
            
         except FileNotFoundError as e:
            # Catch the error and print the specific file that failed
            print(f"--- ERROR: File not found ---")
            print(f"Failed to open file: {full_path}")
            print(f"In case: {cs}")
            # Re-raise the exception or handle it as needed
            raise e
         except Exception as e:
            # Catch other potential errors (like NetCDF decoding issues)
            print(f"--- ERROR: General failure ---")
            print(f"Failed to process file: {full_path}")
            print(f"Details: {e}")
            raise e
   
      # After processing all variables for the current case, append the dictionary to the main list
      trmmc.append(case_data)
   '''

   print('Processing MMC transports...')
   trmmc = [{kk: ux.open_mfdataset(CAMGR, os.path.join(CAMDIR % cs, vv[1])) for kk, vv in tprtfiles.items()} for cs in CASES]
   trmmc = [{kk: list(vv.data_vars.values())[-1] for kk, vv in dd.items()} for dd in trmmc]
   trmmc = [ux.UxDataset(data_vars=dd).groupby('time.month').mean('time') for dd in trmmc]
   mmcds = xr.concat([ds.expand_dims(case=[ALIASES[ii]]) for ii, ds in enumerate(trmmc)], dim='case')

   print('Computing eddy as total minus MMC...')
   eddds = totds - mmcds
   circs = ['total', 'mmc', 'eddy']
   masterds = xr.concat([ds.expand_dims(circ=[circs[ii]]) for ii, ds in enumerate([totds, mmcds, eddds])], dim='circ')
   masterds = ux.UxDataset(masterds, uxgrid=ux.open_grid(CAMGR))
   print(masterds)

   print('Aggregating zonally and seasonally...')
   momds = xr.open_mfdataset(MOMH % CASES[1])
   zmds = masterds.map(lambda uxda: uxda.zonal_mean(lat=momds['yh'].data))
   szds = zmds.assign_coords(coords=dict(month=masterds.month, circ=masterds.circ, case=masterds.case)).map(lambda da: monthly2sznl(da))
   print(szds)
   mirds = szds.map(lambda da: stack_hemi_sznl(da, antisym=True))

   coslat = np.cos(np.deg2rad(mirds.latitudes))
   zontot_transport = mirds * 2 * np.pi * a_e * coslat / g
   zontot_transport.to_netcdf(FILO)

   print(sys.argv[0], 'done.')

def main_plot():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   ds = xr.open_dataset(FILO)
   coslat = np.cos(np.deg2rad(ds.latitudes))
   sinlat = np.sin(np.deg2rad(ds.latitudes))
   ds = ds.assign(variables=dict(LE=lv*ds.q, SE=cp*ds.T, GP=g*ds.Z, AMu=a_e*coslat*ds.u))
   ds = ds.assign(variables=dict(menth=ds.LE+ds.SE, MSE=ds.LE+ds.SE+ds.GP))
   
   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (30, 6)
   subplot_kw = dict(xlim=(-1, 1), sharey=False)
   lncolors = ['blue', 'orange']
   linkw = dict(total=dict(linewidth=2.0, linestyle='solid'), mmc=dict(linewidth=0.8, linestyle='solid'), eddy=dict(linewidth=1.2, linestyle='dotted'))

   for dv in ds.data_vars:
      fig, axes = plt.subplots(1, 3, subplot_kw=subplot_kw)
      fig.suptitle('%s meridional transport' % str(dv))
      for ii, cs in enumerate(ds['case']):
         ax = axes[ii]
         pltda = ds[dv].sel(case=cs)
         if not cs == 'CTRL':
            pltda = pltda - ds[dv].sel(case='CTRL')
         for tt, szn in enumerate(ds['season']):
            for jj, cr in enumerate(ds['circ']):
               ax.plot(sinlat, pltda.sel(season=szn, circ=cr), color=lncolors[tt], **linkw[str(cr.data)])
         ax.axhline(0, c='gray', lw=0.5)
         [ax.axvline(np.sin(np.deg2rad(ll)), c='gray', lw=0.5) for ll in LATLAB]
         ax.set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))

         ylims = ax.get_ylim()
         maxy = max(np.abs(ylims))
         ax.set_ylim(-maxy, maxy)
         if cs == 'SEED':
            ax.set_ylim(*axes[1].get_ylim())

      plt.savefig(os.path.join(DIRO, '%s.png' % str(dv)), bbox_inches='tight')
      plt.close()

#too slow
def merid_transport_no_stationary(vwnd, sclr, vxsclr, dp, lats=(-90, 90, 5)):
   prodofmeans = vwnd * sclr
   fld2tr = lambda fld: 2 * np.pi * a_e * ((fld * dp).sum(dim='lev') / g).zonal_mean(lat=lats) #3D field to vertically integrated transport
   meantr = fld2tr(prodofmeans)
   tottr = None if vxsclr is None else fld2tr(vxsclr)
   return tottr, meantr

if __name__ == '__main__':
   if MODE == 'compute':
      main()
   if MODE == 'plot':
      main_plot()
