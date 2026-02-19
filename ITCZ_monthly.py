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
OLAT = 40.

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YINV = lambda mu: np.rad2deg(np.arcsin(mu))

def main():
   dss = [xr.open_mfdataset(os.path.join(dr, '*.h0a.0005-08*.nc')) for dr in DIRS]
   dss = [ds.assign_coords(coords=dict(mu=YSCL(ds['lat']))) for ds in dss]
   for ii, ds in enumerate(dss):
      if ds[pres_name].units == 'Pa':
         dss[ii] = ds.assign_coords(coords={pres_name: ds[pres_name] / 100})

   psis = [comppsi(ds)[0] for ds in dss]
   seldom = [da.sel(lat=slice(0, OLAT)).interp(plev=mylev) for da in psis]

   pfit = [da.polyfit('mu', 12) for da in seldom]
   #print(pfit[0])

   plt.rcParams['figure.figsize'] = (5, 10)
   fig, axes = plt.subplots(3, 1)
   for ii, ax in enumerate(axes):
      print(xr.polyval(seldom[ii]['mu'], pfit[ii]))
      ax.scatter(seldom[ii]['lat'], seldom[ii], s=0.5)
      ax.plot(seldom[ii]['lat'], xr.polyval(seldom[ii]['mu'], pfit[ii])['polyfit_coefficients'].squeeze())

   plt.show()

if __name__ == '__main__':
   main()
