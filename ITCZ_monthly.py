import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from TEM_sznl import comppsi

DIRS = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm/hist_regrid_0.25x0.25_onpres/',\
        '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_regrid_0.25x0.25_onpres/',\
        '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist_regrid_0.25x0.25_onpres/']
pres_name = 'plev'
mylev = 700.
ILAT, OLAT = 5., 30.
WINDOW = 2. #deglat
STEP = 0.01

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YINV = lambda mu: np.rad2deg(np.arcsin(mu))

def main():
   dss = [xr.open_mfdataset(os.path.join(dr, '*.h0a.000[5-9]*.nc')) for dr in DIRS]
   dss = [ds.sel(time=ds['time'].dt.month.isin(np.arange(6, 12))) for ds in dss]
   dss = [ds.assign_coords(coords=dict(mu=YSCL(ds['lat']))) for ds in dss]
   for ii, ds in enumerate(dss):
      if ds[pres_name].units == 'Pa':
         dss[ii] = ds.assign_coords(coords={pres_name: ds[pres_name] / 100})

   psis = [comppsi(ds)[0].sel(lat=slice(ILAT, OLAT)).interp(plev=mylev) for ds in dss]
   root_init = [np.abs(ps).idxmin(dim='lat') for ps in psis]
   max_init = [ps.idxmax(dim='lat') for ps in psis]

   rt_intp = [ps.interp(lat=np.arange(root_init[ii] - WINDOW, root_init[ii] + WINDOW + STEP, STEP), method='cubic') for ii, ps in enumerate(psis)]
   mx_intp = [ps.interp(lat=np.arange(max_init[ii] - WINDOW, max_init[ii] + WINDOW + STEP, STEP), method='cubic') for ii, ps in enumerate(psis)]
   #mx_intp

   #plt.rcParams['figure.figsize'] = (5, 10)
   #fig, axes = plt.subplots(3, 1, sharex=True, sharey=True)
   #for ii, ax in enumerate(axes):
   #   ax.scatter(psis[ii]['lat'], psis[ii].squeeze(), s=0.5)
   #   ax.plot(rt_intp[ii]['lat'], rt_intp[ii].squeeze(), color='orange')
   #   ax.plot(mx_intp[ii]['lat'], mx_intp[ii].squeeze(), color='green')

   #plt.show()

   rt_fin = [np.abs(ri).idxmin(dim='lat') for ri in rt_intp]
   mx_fin = [mi.idxmax(dim='lat') for mi in mx_intp]
   print(rt_fin)
   print(mx_fin)

   itcz_loc = rt_fin
   itcz_wid = [mx_fin[ii] - rf for ii, rf in enumerate(rt_fin)]

   plt.rcParams['figure.figsize'] = (12, 5)
   fig, axes = plt.subplots(1, 3, sharex=True, sharey=True)
   for ii, ax in enumerate(axes):
      ax.hist(itcz_loc[ii])
   plt.savefig('itcz_loc_test.png', bbox_inches='tight')
   plt.close()

   plt.rcParams['figure.figsize'] = (12, 5)
   fig, axes = plt.subplots(1, 3, sharex=True, sharey=True)
   for ii, ax in enumerate(axes):
      ax.hist(itcz_wid[ii])
   plt.savefig('itcz_wid_test.png', bbox_inches='tight')
   plt.close()

if __name__ == '__main__':
   main()
