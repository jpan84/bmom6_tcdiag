import sys
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from sznl_funcs import monthly2sznl, stack_hemi_sznl

CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
YMM = '/glade/campaign/univ/upsu0032/jpan_aquaptc/%s/atm/hist_regrid_0.25x0.25_onpres/ymonmean.nc'

pltvars = ['CLDICE', 'CLDLIQ', 'Q', 'QRL', 'QRS']
clevs = dict(cldfrac=2.**np.arange(-1, 7, 2) / 100, tdot=0.5 * np.arange(1, 10))
clevs = {k: np.concatenate((-clevs[k][::-1], clevs[k])) for k in clevs}
levfac = [0.1, 1, 0.25]

clevs = dict(cldfrac=np.array([.1, 1, 10, 25]) / 100, tdot=np.array([.01, .1, 1, 2.5]), cldamt=10. ** np.arange(-7, -1))
clevs = {k: np.concatenate((-clevs[k][::-1], clevs[k])) for k in clevs}
qlevs = [np.arange(-3e-4, 3.1e-4, 5e-5), np.arange(0, 3.1e-2, 5e-3), np.arange(-3.5e-3, 3.6e-3, 5e-4)]

def main():
   dss = [xr.open_dataset(YMM % cs).rename(time='month') for cs in CASES]

   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (30, 12)
   latlim = 1 / np.sqrt(2)
   subplot_kw = dict(xlim=(-latlim, latlim), ylim=(200, 1000), yscale='log')
   fig, axes = plt.subplots(2, 3, layout='constrained', sharex=True, sharey=True, subplot_kw=subplot_kw)
   axes[0][0].invert_yaxis()

   sznl_stacked = dict()
   for dv in pltvars:
      sznl = [monthly2sznl(ds[dv].mean(dim='lon')) for ds in dss]
      sznl_stacked[dv] = [stack_hemi_sznl(sz, antisym=False, latnm='lat') for sz in sznl]
   sznl_stacked['CLOUD'] = sznl_stacked['CLDLIQ'] + sznl_stacked['CLDICE']

   for csi in range(len(CASES)):
      for szj, szn in enumerate(['JJA', 'SON']):
         ax, do_diff = axes[szj][csi], csi != 1
         sinlat = np.sin(np.deg2rad(dss[csi]['lat']))
         pltq = sznl_stacked['Q'][csi].sel(season=szn) if not do_diff else sznl_stacked['Q'][csi].sel(season=szn) - sznl_stacked['Q'][1].sel(season=szn)
         csf = ax.contourf(sinlat, dss[csi]['plev'] / 100, pltq, cmap='BrBG' if do_diff else 'YlGnBu', levels=qlevs[csi])
         plt.colorbar(csf, ax=ax)
         pltcld = sznl_stacked['CLOUD'][csi].sel(season=szn) if not do_diff else sznl_stacked['CLOUD'][csi].sel(season=szn) - sznl_stacked['CLOUD'][1].sel(season=szn)
         pltlw = sznl_stacked['QRL'][csi].sel(season=szn) if not do_diff else sznl_stacked['QRL'][csi].sel(season=szn) - sznl_stacked['QRL'][1].sel(season=szn)
         pltsw = sznl_stacked['QRS'][csi].sel(season=szn) if not do_diff else sznl_stacked['QRS'][csi].sel(season=szn) - sznl_stacked['QRS'][1].sel(season=szn)
         ax.contour(sinlat, dss[csi]['plev'] / 100, pltcld, colors='gray', levels=clevs['cldamt'])# * levfac[csi])
         ax.contour(sinlat, dss[csi]['plev'] / 100, 0 * pltlw * 86400, colors='green', levels=clevs['tdot'])# * levfac[csi])
         ax.contour(sinlat, dss[csi]['plev'] / 100, 0 * pltsw * 86400, colors='orange', levels=clevs['tdot'])# * levfac[csi])
         ax.set_xticks(np.sin(np.deg2rad(np.arange(-45, 46, 5))), np.arange(-45, 46, 5))
         ax.set_yticks(np.arange(200, 1001, 100), np.arange(200, 1001, 100))
   plt.savefig('q_cldamt.png', bbox_inches='tight')
   plt.close()

print(sys.argv[0], 'done.')

if __name__ == '__main__':
   main()
