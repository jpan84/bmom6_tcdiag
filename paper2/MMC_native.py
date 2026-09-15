import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post')
from paths import ARCHRT, ALIA, CTLIX, CASENAMES, CAMGR
import consts as c
import uxarray as ux
import xarray as xr
import numpy as np
from sznl_funcs import stack_hemi_sznl, monthly2sznl, agg_time
import matplotlib.pyplot as plt

h0 = '/glade/campaign/univ/upsu0032/jpan_aquaptc/%s/atm/hist/*.h0a.0009-1*.nc'

def main_compute():
   dss = [ux.open_mfdataset(CAMGR, h0 % cs).expand_dims(case=[ALIA[ii]]) for ii, cs in enumerate(CASENAMES)]
   for ii in range(1, len(dss)):
      dss[ii].uxgrid = dss[0].uxgrid
   ds = ux.concat(dss, dim='case')

   print('Setting up coords...')
   aterm = ds['hyai'] * c.P0
   bterm = ds['hybi'] * ds['PS']
   p_ilev = aterm + bterm
   dp3d = p_ilev.diff('ilev').rename(dict(ilev='lev')).assign_coords(lev=ds['lev'])
   ds = ds.assign(variables=dict(dp3d=dp3d, p_ilev=p_ilev, coslat=np.cos(np.deg2rad(ds['lat']))))

   vdp = ds['V'] * ds['dp3d']
   vdp_zm = vdp.zonal_mean((-90, 90, 1.5)).assign_coords(case=ds['case'], time=ds['time'])
   p_ilev_zm = ds['p_ilev'].zonal_mean((-90, 90, 1.5)) #the integrals values correspond to interfaces

   #cum = vdp_zm.cumulative('lev').sum()
   vdp_zm_ilev = vdp_zm.rename(lev='ilev')
   
   vdp_zm_ilev = xr.concat(
       [xr.zeros_like(vdp_zm_ilev.isel(ilev=0)), vdp_zm_ilev],
       dim='ilev',
   )
   
   vdp_zm_ilev = vdp_zm_ilev.assign_coords(ilev=p_ilev_zm.ilev)
   
   cum = vdp_zm_ilev.cumsum('ilev')
   p_tgt = np.concatenate((ds['lev'].data, [1010.])) * 100.

   myint = lambda data, pi: np.interp(p_tgt, pi, data, left=np.nan, right=np.nan)

   cum = vdp_zm_ilev.cumsum('ilev').chunk({'ilev': -1})
   p_ilev_zm = p_ilev_zm.chunk({'ilev': -1})
   print(p_ilev_zm.isel(ilev=-1).max().values, p_ilev_zm.isel(ilev=-1).min().values)

   cum_p = xr.apply_ufunc(
           myint, cum, p_ilev_zm,
           input_core_dims=[['ilev'], ['ilev']],
           output_core_dims=[['plev']],
           vectorize=True,
           dask='parallelized',
           output_dtypes=[cum.dtype],
           dask_gufunc_kwargs={'output_sizes': {'plev': len(p_tgt)}},
   ).assign_coords(plev=p_tgt)

   print(cum_p)

   cum_p.to_dataset.to_netcdf('test_MMC_native.nc')

   exit()


if __name__ == '__main__':
   main_compute()
