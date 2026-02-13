import sys
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from sznl_funcs import monthly2sznl, stack_hemi_sznl

CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
YMM = '/glade/campaign/univ/upsu0032/jpan_aquaptc/%s/atm/hist_regrid_0.25x0.25_onpres/ymonmean.nc'

letters = [['(a)', '(b)', '(c)'], ['(d)', '(e)', '(f)']]
CSTTL = ['UNSEED$-$CTRL', 'CTRL', 'SEED$-$CTRL']

pltvars = ['T', 'THETA']
temp_levs = [np.arange(-0.5, 0.51, 0.05), np.arange(210, 316, 5), np.arange(-2, 2.1, 0.25)]
thta_levs = [temp_levs[0], np.arange(270, 401, 10), temp_levs[-1]]
thta_levs = [tl[tl != 0] for tl in thta_levs]

def main():
   dss = [xr.open_dataset(YMM % cs).rename(time='month') for cs in CASES]
   dss = [ds.assign(THETA=ds['T'] * (1e5 / ds['plev'])**(287/1005)) for ds in dss]

   plt.rc('font', size=20)
   plt.rcParams['figure.figsize'] = (30, 12)
   latlim = np.sin(np.deg2rad(50))#1 / np.sqrt(2)
   subplot_kw = dict(xlim=(-latlim, latlim), ylim=(200, 1000), yscale='log')
   fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, subplot_kw=subplot_kw) #layout='constrained'
   axes[0][0].invert_yaxis()

   sznl_stacked = dict()
   for dv in pltvars:
      sznl = [monthly2sznl(ds[dv].mean(dim='lon')) for ds in dss]
      sznl_stacked[dv] = [stack_hemi_sznl(sz, antisym=False, latnm='lat') for sz in sznl]


   for csi in range(len(CASES)):
      for szj, szn in enumerate(['JJA', 'SON']):
         ax, do_diff = axes[szj][csi], csi != 1
         sinlat = np.sin(np.deg2rad(dss[csi]['lat']))
         pltT = sznl_stacked['T'][csi].sel(season=szn) if not do_diff else sznl_stacked['T'][csi].sel(season=szn) - sznl_stacked['T'][1].sel(season=szn)
         csf = ax.contourf(sinlat, dss[csi]['plev'] / 100, pltT, cmap='bwr' if do_diff else 'YlGnBu', levels=temp_levs[csi])
         plt.colorbar(csf, ax=ax)
         pltth = sznl_stacked['THETA'][csi].sel(season=szn) if not do_diff else sznl_stacked['THETA'][csi].sel(season=szn) - sznl_stacked['THETA'][1].sel(season=szn)
         ax.contour(sinlat, dss[csi]['plev'] / 100, pltth, colors='black', levels=thta_levs[csi])# * levfac[csi])
         ax.contour(sinlat, dss[csi]['plev'] / 100, sznl_stacked['T'][csi].sel(season=szn), levels=[273.15], colors='purple')
         ax.contour(sinlat, dss[csi]['plev'] / 100, sznl_stacked['T'][1].sel(season=szn), levels=[273.15], colors='brown')

         ax.tick_params(right=True, top=True, labelbottom=True, labelleft=True)
         ax.set_xticks(np.sin(np.deg2rad(np.arange(-50, 51, 5))), np.arange(-50, 51, 5))
         [lbl.set_visible(False) for ii, lbl in enumerate(ax.get_xticklabels()) if ii % 2]
         ax.set_yticks(np.arange(200, 1001, 100), np.arange(200, 1001, 100))
         ax.set_title(letters[szj][csi], loc='left')
         ax.set_xlabel('Latitude')

         if szj == 0:
            ax.set_title(CSTTL[csi], fontsize=28)
         if csi == 0:
            ax.set_ylabel(szn + '   ', rotation='horizontal', fontsize=28)
   plt.savefig('T_theta.png', bbox_inches='tight')
   plt.show()
   plt.close()

   print(sys.argv[0], 'done.')

if __name__ == '__main__':
   main()
