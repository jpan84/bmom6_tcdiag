import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import numpy as np

DIRI = './streamf_sznl_seedmatch'
FILI = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl__TEM.nc'

#DIRI = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr/atm'
#FILI = '0012_JJASON_onpres_0.25_EPF_tm.nc'
DIFF = False

FILI = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl_b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch_TEM.nc'
DIFF = True

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-90, 90.1, 30).astype(np.int_)
YLOC = YSCL(YLAB)

def main():
   ds = xr.open_dataset(os.path.join(DIRI, FILI))

   plotfields = [ds[dv] for dv in ['EPy_EMF_d', 'EPy_adv_d', 'EPz_EHF_d', 'EPz_adv_d']]
   plotfields.append(sum(plotfields))
   plottitles = ['EPFd (y) due to EMF', 'EPFd (y) due to mean wu adv', 'EPFd (z) due to EHF', 'EPFd (z) due to mean vu adv', 'Sum']

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
         curplt = plotfields[sfi] * 86400 # * (1 / 86400 / 100)
         if 'season' in curplt.dims:
            curplt = curplt.isel(season=szj)
         CSF = ax.contourf(YSCL(ds['lat']), ds['plev'], curplt, **contourfkwargs)
         CS1 = ax.contour(YSCL(ds['lat']), ds['plev'], curplt, **contourkwargs)
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

   print('Plotting EPF and EPFd...')
   plt.rcParams['figure.figsize'] = (10, 12)
   fig, axes = plt.subplots(2, 1, layout='constrained', sharex=True, sharey=True, subplot_kw=subplot_kw)
   EPy = ds['EPy_EMF'] + ds['EPy_adv']
   EPz = ds['EPz_EHF'] + ds['EPz_adv']

   islc = dict(lat=slice(None, None, 16), plev=slice(None, -2))
   for szj, szn in enumerate(['JJA', 'SON']):
      ax = axes[szj]
      divplt = plotfields[-1].isel(**islc) * 86400 # * (1 / 86400 / 100)
      pltepy, pltepz = EPy.isel(**islc) / 1e2 / 1e1, EPz.isel(**islc)
      pltlat, pltp = ds['lat'].isel(lat=islc['lat']), ds['plev'].isel(plev=islc['plev'])
      if 'season' in divplt.dims:
         divplt = divplt.isel(season=szj)
         pltepy, pltepz = pltepy.sel(season=szn), pltepz.sel(season=szn)
      csf = ax.contourf(pltlat, pltp, divplt, **contourfkwargs)
      #CS1 = ax.contour(YSCL(ds['lat']), ds['plev'], curplt, **contourkwargs)
      qv = ax.quiver(pltlat, pltp, pltepy, pltepz, scale=1e6, pivot='mid') #scale=5e7
      if szj == 0:
         ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
         ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
         ax.set_ylabel('Pressure [hPa]')
         ax.set_yticks([50, 70, 100, 200, 300, 500, 700, 850], [50, 70, 100, 200, 300, 500, 700, 850])#np.arange(100, 1001, 100))
         cb = plt.colorbar(csf, ax=axes)
         cbt = cb.get_ticks()
         cbt = np.concatenate((cbt[:cbt.size // 2 + 1], -cbt[:cbt.size // 2 + 1][::-1]))
         cb.set_ticks(cbt)
         cb.set_ticklabels(['%.1f' % tk for tk in cbt])
         ax.invert_yaxis()
      ax.set_xticks(np.arange(-90, 91, 10))
      plt.xlim(-10, 60)
      ax.set_title('%s\n[m s$^{-1}$ day$^{-1}$]' % plottitles[sfi])
      ax.set_xlabel('Latitude [°]')

   #fig.tight_layout()
   plt.savefig(os.path.join(DIRI, '%s_EPF_EPFd.png' % FILI.split('.')[-2]), bbox_inches='tight')
   plt.close()

if __name__ == '__main__':
   main()
