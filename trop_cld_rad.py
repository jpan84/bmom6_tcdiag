import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from sznl_funcs import monthly2sznl, stack_hemi_sznl

CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
YMM = '/glade/derecho/scratch/jpan/archive/%s/atm/hist_regrid_0.25x0.25_onpres/ymonmean.nc'

pltvars = ['CLOUD', 'Q', 'QRL', 'QRS']
clevs = dict(cldfrac=2.**np.arange(-1, 7, 2) / 100, tdot=0.5 * np.arange(1, 10))
clevs = {k: np.concatenate((-clevs[k][::-1], clevs[k])) for k in clevs}
levfac = [0.1, 1, 0.25]

def main():
   dss = [xr.open_dataset(YMM % cs).rename(time='month') for cs in CASES]

   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (30, 12)
   subplot_kw = dict(xlim=(-0.5, 0.5), ylim=(100, 1000), yscale='log')
   fig, axes = plt.subplots(2, 3, layout='constrained', sharex=True, sharey=True, subplot_kw=subplot_kw)
   axes[0][0].invert_yaxis()

   sznl_stacked = dict()
   for dv in pltvars:
      sznl = [monthly2sznl(ds[dv].mean(dim='lon')) for ds in dss]
      sznl_stacked[dv] = [stack_hemi_sznl(sz, antisym=False, latnm='lat') for sz in sznl]

   for csi in range(len(CASES)):
      for szj, szn in enumerate(['JJA', 'SON']):
         ax, do_diff = axes[szj][csi], csi != 1
         sinlat = np.sin(np.deg2rad(dss[csi]['lat']))
         pltq = sznl_stacked['Q'][csi].sel(season=szn) if not do_diff else sznl_stacked['Q'][csi].sel(season=szn) - sznl_stacked['Q'][1].sel(season=szn)
         csf = ax.contourf(sinlat, dss[csi]['plev'] / 100, pltq, cmap='BrBG' if do_diff else 'YlGnBu')
         plt.colorbar(csf, ax=ax)
         pltcld = sznl_stacked['CLOUD'][csi].sel(season=szn) if not do_diff else sznl_stacked['CLOUD'][csi].sel(season=szn) - sznl_stacked['CLOUD'][1].sel(season=szn)
         pltlw = sznl_stacked['QRL'][csi].sel(season=szn) if not do_diff else sznl_stacked['QRL'][csi].sel(season=szn) - sznl_stacked['QRL'][1].sel(season=szn)
         pltsw = sznl_stacked['QRS'][csi].sel(season=szn) if not do_diff else sznl_stacked['QRS'][csi].sel(season=szn) - sznl_stacked['QRS'][1].sel(season=szn)
         ax.contour(sinlat, dss[csi]['plev'] / 100, pltcld, colors='gray', levels=clevs['cldfrac'] * levfac[csi])
         ax.contour(sinlat, dss[csi]['plev'] / 100, pltlw * 86400, colors='green', levels=clevs['tdot'] * levfac[csi])
         ax.contour(sinlat, dss[csi]['plev'] / 100, pltsw * 86400, colors='orange', levels=clevs['tdot'] * levfac[csi])
         ax.set_xticks(np.sin(np.deg2rad(np.arange(-30, 31, 5))), np.arange(-30, 31, 5))
   plt.savefig('q_cld_rad.png', bbox_inches='tight')
   plt.close()

if __name__ == '__main__':
   main()
