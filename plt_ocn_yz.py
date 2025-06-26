import xarray as xr
import numpy as np
import os
import sznl_funcs
import matplotlib.pyplot as plt
import matplotlib.colors as colors

HISTS = '/glade/derecho/scratch/jpan/archive/%s/ocn/hist/*mom6.hm*[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'
CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
ALIASES = ['UNSEED', 'CTRL', 'SEED']
DIRO = './ocn_yz_plts/'

DIMS3D = {('time', 'zl', 'yq', 'xq'), ('time', 'zl', 'yq', 'xh'), ('time', 'zl', 'yh', 'xq'), ('time', 'zl', 'yh', 'xh')}
ANTISYM_STR = ['vorticity', 'meridional', 'y transport', 'y velocity']

ZSCL = np.vectorize(lambda zl: zl / 2850 if zl <= 2000 else 2000 / 2850 + (zl - 2000) / 2000 * 850 / 2850)
ZLAB = np.arange(1000, 5000, 1000)
ZLOC = ZSCL(ZLAB)

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-90, 90.1, 30).astype(np.int_)
YLOC = YSCL(YLAB)

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   dss = [xr.open_mfdataset(HISTS % cs) for cs in CASES]
   dvars = list(dss[0].data_vars)

   for dv in dvars:
      if dss[0][dv].dims not in DIMS3D:
         print('%s not 3D. Skipping...' % dv)
         continue

      print('Working on %s...' % dv)
      xmeth = dss[0][dv].attrs['cell_methods'].split()[-3].split(':') #e.g., ['xh', 'mean']
      ymeth = dss[0][dv].attrs['cell_methods'].split()[-4].split(':')
      zmeth = dss[0][dv].attrs['cell_methods'].split()[-5].split(':')

      zonagg = None
      if xmeth[1] in {'point', 'mean'}:
         zonagg = [ds[dv].mean(dim=xmeth[0]) for ds in dss]
      else:
         zonagg = [ds[dv].sum(dim=xmeth[0]) for ds in dss]

      if zmeth[1] == 'sum':
         zonagg = [za.cumsum(dim=zmeth[0]) for za in zonagg]

      monmeans = [za.groupby('time.month').mean() for za in zonagg]
      if ymeth[0] == 'yq':
         monmeans = [mm.isel(yq=slice(None, -1)) for mm in monmeans]
      antisym = False
      for sstr in ANTISYM_STR:
         if sstr in dss[0][dv].attrs['long_name'].lower():
            antisym = True; break
      sznl = [sznl_funcs.stack_hemi_sznl(sznl_funcs.monthly2sznl(mm), antisym=antisym, latnm=ymeth[0]) for mm in monmeans]
      vmax = max([np.abs(sz).max().values for sz in sznl])
      vmin = min([sz.min().values for sz in sznl])
      zero_centered = True
      if vmax - vmin < 0.25 * vmax:
         zero_centered = False

      plt.rc('font', size=16)
      plt.rcParams['figure.figsize'] = (30, 12)
      contourfkwargs = {'cmap': 'coolwarm' if zero_centered else 'rainbow', 'norm': colors.CenteredNorm(vcenter=0, halfrange=vmax) if zero_centered else None}
      subplot_kw = dict(xlim=(-1, 1), ylim=(0, 1))
      fig, axes = plt.subplots(2, 3, layout='constrained', sharey=True, subplot_kw=subplot_kw)
      plt.suptitle(f"{dv} [{dss[0][dv].attrs['units']}]\n{dss[0][dv].attrs['long_name']}\n{dss[0][dv].attrs['cell_methods']}")

      for ii, sz in enumerate(['JJA', 'SON']):
         for jj, cs in enumerate(CASES):
            ax = axes[ii, jj]
            csf = ax.contourf(YSCL(sznl[jj][ymeth[0]]), ZSCL(sznl[jj][zmeth[0]]), sznl[jj].sel(season=sz).data, **contourfkwargs)
            if ii == 0 and jj == 0:
               ax.set_ylabel('Depth [m]')
               ax.set_yticks(ZLOC, ZLAB)
               ax.invert_yaxis()
               plt.colorbar(csf, ax=axes)
            ax.set_xticks(YLOC, YLAB)
            ax.set_xlabel('Latitude [°]')

      plt.savefig(os.path.join(DIRO, '%s.png' % dv), bbox_inches='tight')
      plt.close()


if __name__ == '__main__':
   main()
