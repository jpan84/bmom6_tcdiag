import sys
import os
import math
import numpy as np
import xarray as xr
import uxarray as ux
from dask.diagnostics import ProgressBar
from sznl_funcs import monthly2sznl, stack_hemi_sznl

MODE = sys.argv[1] #'compute', 'plot'

FILO = 'transports_3D_conserv.nc'

CAMGR = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
CAMDIR = '/glade/campaign/univ/upsu0032/jpan_aquaptc//%s/atm/hist/'
MOMH = '/glade/derecho/scratch/jpan/archive/%s/ocn/hist/*mom6.hm*[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'
CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
ALIASES = ['UNSEED', 'CTRL', 'SEED']


P0 = 1e5
g = 9.79764
OM = 7.2921e-5
a_e = 6.371e6

#inflds = [('Q', 'VQ')]#, ('T', 'VT'), ('U', 'VU'), ('Z3', None)]
tprtfiles = dict(q=('vqdp_monthly.nc', 'vqdp_mmc_monthly.nc'), T=('vTdp_monthly.nc', 'vTdp_mmc_monthly.nc'), u=('vudp_monthly.nc', 'vudp_mmc_monthly.nc'), Z=(None, 'vZdp_mmc_monthly.nc'))

def main():
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

   trmmc = [{kk: ux.open_mfdataset(CAMGR, os.path.join(CAMDIR % cs, vv[1])) for kk, vv in tprtfiles.items()} for cs in CASES]
   trmmc = [{kk: list(vv.data_vars.values())[-1] for kk, vv in dd.items()} for dd in trmmc]
   trmmc = [ux.UxDataset(data_vars=dd).groupby('time.month').mean('time') for dd in trmmc]
   mmcds = xr.concat([ds.expand_dims(case=[ALIASES[ii]]) for ii, ds in enumerate(trmmc)], dim='case')

   eddds = totds - mmcds
   circs = ['total', 'mmc', 'eddy']
   dss = xr.concat([ds.expand_dims(circ=circs[ii]) for ii, ds in enumerate([totds, mmcds, eddds])], dim='circ')
   print(dss)


def main_slow():
   print('Opening history files...')
   camdss = [ux.open_mfdataset(CAMGR, CAMH % cs) for cs in CASES]
   momdss = [xr.open_mfdataset(MOMH % cs) for cs in CASES]

   print('Obtaining MOM lats and dp3d...')
   mom_h_lat = momdss[0]['yh']
   aterm = (camdss[0]['hyai'].isel(ilev=slice(1, None)) - camdss[0]['hyai'].isel(ilev=slice(None, -1)).data) * P0
   bterm = (camdss[0]['hybi'].isel(ilev=slice(1, None)) - camdss[0]['hybi'].isel(ilev=slice(None, -1)).data) * camdss[0]['PS']
   dp3d = aterm + bterm
   dp3d = dp3d.rename(dict(ilev='lev')).assign_coords(lev=camdss[0]['lev'])

   #print('Computing V Q transport...')
   #vqtr = [merid_transport_no_stationary(ds['V'], ds['Q'], ds['VQ'], dp3d, lats=mom_h_lat.data) for ds in camdss]
   #vqdss = [xr.Dataset(data_vars=dict(vqmean=tup[1].to_xarray(), vqtot=tup[0].to_xarray())) for tup in vqtr]
   #print(vqdss[0])

   outdss = [xr.Dataset(coords=dict(case=ali)) for ali in ALIASES]
   for ii, fldtup in enumerate(inflds):
      print('Working on transports', fldtup, '...')
      tprt = [merid_transport_no_stationary(ds['V'], ds[fldtup[0]], None if fldtup[1] is None else ds[fldtup[1]], dp3d, lats=mom_h_lat.data) for ds in camdss]
      tprtdss = [xr.Dataset(data_vars={'V' + fldtup[0] + '_mean': tup[1].to_xarray(), 'V' + fldtup[0] + '_tot': tup[0].to_xarray()}).assign_coords(\
                   coords=dict(time=camdss[jj]['time'])) for jj, tup in enumerate(tprt)]
      momeans = [td.groupby('time.month').mean('time') for td in tprtdss]
      szmeans = [mm.map(lambda da: stack_hemi_sznl(monthly2sznl(da), antisym=True)) for mm in momeans]
      outdss = [xr.merge([ds, szmeans[jj]]) for jj, ds in enumerate(outdss)]

   outds = xr.concat(outdss, dim='case')
   print('Saving output .nc...')
   with ProgressBar():
      outds.to_netcdf(FILO)
   print(sys.argv[0], 'done.')

def main_plot():
   return

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
