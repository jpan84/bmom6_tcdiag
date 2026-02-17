import sys
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from sznl_funcs import monthly2sznl, stack_hemi_sznl

cp = 1005
Lv = 2.501e6
g = 9.807

YMM = '/glade/campaign/univ/upsu0032/jpan_aquaptc/%s/atm/hist_regrid_0.25x0.25_onpres/ymonmean.nc'
YMM2 = '/glade/derecho/scratch/jpan/archive/%s/atm/hist_regrid_0.25x0.25_onpres/ymonmean.nc'
CASES = [(YMM, 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed'), (YMM, 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl'), (YMM2, 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch')]

letters = [['(a)', '(b)', '(c)'], ['(d)', '(e)', '(f)']]
CSTTL = ['UNSEED$-$CTRL', 'CTRL', 'MSEED$-$CTRL']

pltvars = ['SE', 'LE', 'GP', 'DSE', 'MSE']
clevs = dict(MSE=(np.arange(-1000, 1001, 100), np.arange(3e5, 3.9e5, 5e3), np.arange(-1000, 1001, 100)))
clevs['DSE'] = clevs['MSE']


temp_levs = [np.arange(-0.5, 0.51, 0.05), np.arange(210, 316, 5), np.arange(-2, 2.1, 0.25)]
temp_levs = [np.arange(-0.5, 0.51, 0.05), np.arange(210, 316, 5), np.arange(-0.5, 0.51, 0.05)]
thta_levs = [temp_levs[0], np.arange(270, 401, 10), temp_levs[-1]]
thta_levs = [tl[tl != 0] for tl in thta_levs]

def main():
   dss = [xr.open_dataset(cs[0] % cs[1]).rename(time='month') for cs in CASES]
   dss = [ds.assign(SE=cp * ds['T']) for ds in dss]
   dss = [ds.assign(LE=Lv * ds['Q']) for ds in dss]
   dss = [ds.assign(GP=g * ds['Z3']) for ds in dss]
   dss = [ds.assign(DSE=ds['SE'] + ds['GP']) for ds in dss]
   dss = [ds.assign(MSE=ds['DSE'] + ds['LE']) for ds in dss]

   plt.rc('font', size=20)
   plt.rcParams['figure.figsize'] = (30, 12)
   latlim = np.sin(np.deg2rad(50))#1 / np.sqrt(2)
   subplot_kw = dict(xlim=(-latlim, latlim), ylim=(100, 1000), yscale='log')
   fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, subplot_kw=subplot_kw) #layout='constrained'
   axes[0][0].invert_yaxis()

   sznl_stacked = dict()
   for dv in pltvars:
      sznl = [monthly2sznl(ds[dv].mean(dim='lon')) for ds in dss]
      sznl_stacked[dv] = [stack_hemi_sznl(sz, antisym=False, latnm='lat') for sz in sznl]


   myvar = 'DSE'
   for csi in range(len(CASES)):
      for szj, szn in enumerate(['JJA', 'SON']):
         ax, do_diff = axes[szj][csi], csi != 1
         sinlat = np.sin(np.deg2rad(dss[csi]['lat']))

         #pltT = sznl_stacked['T'][csi].sel(season=szn) if not do_diff else sznl_stacked['T'][csi].sel(season=szn) - sznl_stacked['T'][1].sel(season=szn)
         #csf = ax.contourf(sinlat, dss[csi]['plev'] / 100, pltT, cmap='bwr' if do_diff else 'YlGnBu', levels=temp_levs[csi])
         #plt.colorbar(csf, ax=ax)
         #pltth = sznl_stacked['THETA'][csi].sel(season=szn) if not do_diff else sznl_stacked['THETA'][csi].sel(season=szn) - sznl_stacked['THETA'][1].sel(season=szn)
         #ax.contour(sinlat, dss[csi]['plev'] / 100, pltth, colors='black', levels=thta_levs[csi])# * levfac[csi])
         #ax.contour(sinlat, dss[csi]['plev'] / 100, sznl_stacked['T'][csi].sel(season=szn), levels=[273.15], colors='purple')
         #ax.contour(sinlat, dss[csi]['plev'] / 100, sznl_stacked['T'][1].sel(season=szn), levels=[273.15], colors='brown')

         pltfld = sznl_stacked[myvar][csi].sel(season=szn) if not do_diff else sznl_stacked[myvar][csi].sel(season=szn) - sznl_stacked[myvar][1].sel(season=szn)
         csf = ax.contourf(sinlat, dss[csi]['plev'] / 100, pltfld, cmap='bwr' if do_diff else 'YlGnBu', levels=clevs[myvar][csi])
         plt.colorbar(csf, ax=ax)

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
   plt.savefig('DSE_test.png', bbox_inches='tight')
   plt.show()
   plt.close()

   print(sys.argv[0], 'done.')

if __name__ == '__main__':
   main()
