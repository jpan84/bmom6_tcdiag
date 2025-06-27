import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import numpy as np

DIRI = './streamf_sznl'
FILI = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl__TEM.nc'
DIFF = False

FILI = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl_b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed_TEM.nc'
DIFF = True

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-90, 90.1, 30).astype(np.int_)
YLOC = YSCL(YLAB)

def main():
   ds = xr.open_dataset(os.path.join(DIRI, FILI))

   plotfields = [ds[dv] for dv in ['EPy_EMF_d', 'EPy_adv_d', 'EPz_EHF_d', 'EPz_adv_d']]
   plotfields.append(sum(plotfields))
   plottitles = ['EPFd (y) due to EMF', 'EPFd (y) due to mean vu adv', 'EPFd (z) due to EHF', 'EPFd (z) due to mean wu adv', 'Sum']

   print('Plotting EPFd...')
   plt.rc('font', size=20)
   plt.rcParams['figure.figsize'] = (52, 12)
   contourkwargs = {'colors': 'black', 'levels': .2 * 2.**np.arange(-1, 7, 1)}
   contourkwargs['levels'] = np.concatenate((-contourkwargs['levels'][::-1], contourkwargs['levels']))
   contourfkwargs = {'cmap': 'RdYlBu_r' if DIFF else 'RdBu_r', 'levels': contourkwargs['levels'], 'norm': colors.SymLogNorm(0.1), 'extend': 'both'} #Use RdYlBu_r for diff
   #clabelkwargs = {'inline': 1, 'fontsize': 10, 'colors': 'black', 'fmt': '%.1f'}
   subplot_kw = dict(xlim=(-1, 1), ylim=(100, 1000), yscale='log')
   fig, axes = plt.subplots(2, 5, layout='constrained', sharey=True, subplot_kw=subplot_kw)

   for sfi, sfn in enumerate(plotfields):
      for szj, szn in enumerate(['JJA', 'SON']):
         ax = axes[szj, sfi]
         CSF = ax.contourf(YSCL(ds['lat']), ds['plev'], plotfields[sfi].isel(season=szj).data * 86400, **contourfkwargs)
         CS1 = ax.contour(YSCL(ds['lat']), ds['plev'], plotfields[sfi].isel(season=szj).data * 86400, **contourkwargs)
         if sfi == 0 and szj == 0:
            ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
            ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
            ax.set_ylabel('Pressure [hPa]')
            ax.set_yticks(np.arange(100, 1001, 100))
            cb = plt.colorbar(CSF, ax=axes)
            cbt = cb.get_ticks()
            cbt = np.concatenate((cbt[:cbt.size // 2 + 1], -cbt[:cbt.size // 2 + 1][::-1]))
            cb.set_ticks(cbt)
            cb.set_ticklabels(['%.1f' % tk for tk in cbt])
            ax.invert_yaxis()
         ax.set_xticks(YLOC, YLAB)
         ax.set_title('%s\n[m s$^{-1}$ day$^{-1}$]' % plottitles[sfi])
         ax.set_xlabel('Latitude [°]')

   #fig.tight_layout()
   plt.savefig(os.path.join(DIRI, '%s_EPFd_terms.png' % FILI.split('.')[-2]), bbox_inches='tight')
   plt.close()

if __name__ == '__main__':
   main()
