import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from TEM_sznl import comppsi

DIRS = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm/hist_regrid_0.25x0.25_onpres/',\
        '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_regrid_0.25x0.25_onpres/',\
        '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist_regrid_0.25x0.25_onpres/']
ALIS = ['UNSEED', 'CTRL', 'MSEED']
pres_name = 'plev'
mylev = 700.
ILAT, OLAT = 5., 30.
WINDOW = 1. #deglat
STEP = 0.02

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YINV = lambda mu: np.rad2deg(np.arcsin(mu))

def main():
   dss = [xr.open_mfdataset(os.path.join(dr, '*.h0a.000[5-9]*.nc')) for dr in DIRS]
   dss = [ds.sel(time=ds['time'].dt.month.isin(np.arange(6, 12))) for ds in dss]
   #dss = [ds.assign_coords(coords=dict(mu=YSCL(ds['lat']))) for ds in dss]
   for ii, ds in enumerate(dss):
      if ds[pres_name].units == 'Pa':
         dss[ii] = ds.assign_coords(coords={pres_name: ds[pres_name] / 100})

   psis = [comppsi(ds)[0].sel(lat=slice(ILAT, OLAT)).interp(plev=mylev) for ds in dss]
   root_init = [np.abs(ps).idxmin(dim='lat') for ps in psis]
   max_init = [ps.idxmax(dim='lat') for ps in psis]

   rt_intp = [ps.interp(lat=np.arange(root_init[ii].min(dim='time') - WINDOW, root_init[ii].max(dim='time') + WINDOW + STEP, STEP), method='cubic') for ii, ps in enumerate(psis)]
   mx_intp = [ps.interp(lat=np.arange(max_init[ii].min(dim='time') - WINDOW, max_init[ii].max(dim='time') + WINDOW + STEP, STEP), method='cubic') for ii, ps in enumerate(psis)]
   #mx_intp

   #plt.rcParams['figure.figsize'] = (5, 10)
   #fig, axes = plt.subplots(3, 1, sharex=True, sharey=True)
   #for ii, ax in enumerate(axes):
   #   ax.scatter(psis[ii]['lat'], psis[ii].squeeze(), s=0.5)
   #   ax.plot(rt_intp[ii]['lat'], rt_intp[ii].squeeze(), color='orange')
   #   ax.plot(mx_intp[ii]['lat'], mx_intp[ii].squeeze(), color='green')

   #plt.show()

   rt_fin = [np.abs(ri).idxmin(dim='lat').expand_dims(case=[ALIS[ii]]) for ii, ri in enumerate(rt_intp)]
   mx_fin = [mi.idxmax(dim='lat').expand_dims(case=[ALIS[ii]]) for ii, mi in enumerate(mx_intp)]
   rt_da = xr.concat(rt_fin, dim='case')
   mx_da = xr.concat(mx_fin, dim='case')
   print(rt_da)
   print(mx_da)

   itcz_loc = rt_da
   itcz_wid = mx_da - rt_da #[mx_fin[ii] - rf for ii, rf in enumerate(rt_fin)]
   outds = xr.Dataset(data_vars=dict(zerolat=rt_da, maxlat=mx_da, itcz_width=itcz_wid))
   outds.to_netcdf('ITCZ_monthly_0005-0009.nc')

   plt.rcParams['figure.figsize'] = (12, 5)
   fig, axes = plt.subplots(1, 3, sharex=True, sharey=True)
   for ii, ax in enumerate(axes):
      ax.hist(itcz_loc.isel(case=ii))
   plt.savefig('itcz_loc_test.png', bbox_inches='tight')
   plt.close()

   plt.rcParams['figure.figsize'] = (12, 5)
   fig, axes = plt.subplots(1, 3, sharex=True, sharey=True)
   for ii, ax in enumerate(axes):
      ax.hist(itcz_wid.isel(case=ii))
   plt.savefig('itcz_wid_test.png', bbox_inches='tight')
   plt.close()

if __name__ == '__main__':
   main()
