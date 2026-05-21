import sys
import xarray as xr
from dask.diagnostics import ProgressBar
import numpy as np
import os
import sznl_funcs
import matplotlib.pyplot as plt
import matplotlib.colors as colors

mode = sys.argv[1]
MODES = ['compute', 'plot']

DO_DIFF = True
HISTS = 'ocn/hist/*mom6.hm*[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'
CASES = ['/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250702_unseed2hPa6m/', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/', '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/']
ALIASES = ['UNSEED_EX', 'UNSEED', 'CTRL', 'SEED', 'SEED_EX']
DIRO = './ocn_yz_plts_5exp_inprog/'

DIMS3D = {('time', 'zl', 'yq', 'xq'), ('time', 'zl', 'yq', 'xh'), ('time', 'zl', 'yh', 'xq'), ('time', 'zl', 'yh', 'xh')}
ANTISYM_STR = ['vorticity', 'meridional', 'y transport', 'y velocity']
INCL = {'thetao', 'vmo', 'vhGM'}

ZSCL = np.vectorize(lambda zl: zl / 2850 if zl <= 2000 else 2000 / 2850 + (zl - 2000) / 2000 * 850 / 2850)
ZLAB = np.arange(1000, 5000, 1000)
ZLOC = ZSCL(ZLAB)

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)

CLEVS = dict(thetao=((0, 37, 4), (-0.5, 0.51, .05)), rhopot0=((1018, 1030, 1), (-.2, .21, .02)), so=((32, 38, .5), (-.2, .2, .02)), T_ady=((-1.6e16, 1.7e16, 4e15), (-2e15, 2.1e15, 2e14)),\
             uo=((-2.25, 0.6, .25), (-.1, .21, .02)))
ZEROCTR = {'vmo', 'uo', 'vo'}
LIMDEPTH = {'uo', 'thetao', 'rhopot0', 'so'}
LIMLAB = np.arange(200, 1500, 200)
LIMLOC = ZSCL(LIMLAB)

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   dss = [xr.open_mfdataset(os.path.join(cs, HISTS)) for cs in CASES]
   #dss = [ds.assign(variables=dict(vmo_resid=\
   #       (ds['vmo'] + ds['vhGM']).assign_attrs(long_name='Ocean Mass Residual Y Transport',\
   #       cell_methods=ds['vmo'].attrs['cell_methods'], units=ds['vmo'].attrs['units'])\
   #       )) for ds in dss]
   dvars = list(dss[0].data_vars)

   outds = None
   for dv in dvars:
      #!COLORLEVS
      #if dv not in ['vmo_resid']:
      #   continue

      if dss[0][dv].dims not in DIMS3D:
         print('%s not 3D. Skipping...' % dv)
         continue

      if INCL is not None and str(dv) not in INCL:
         print('%s excluded. Skipping...' % dv)
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
      with ProgressBar():
         sznl = [sz.compute() for sz in sznl]

      outda = xr.concat(sznl, dim=xr.DataArray(ALIASES, dims='case'))
      if outds is None:
         outds = xr.Dataset(data_vars={dv: outda})
      else:
         outds = outds.assign(variables={dv: outda})

      vmin, vmax, zero_centered = None, None, False
      if DO_DIFF:
         zero_centered = True
         sznl = [(sznl[ii] - sznl[1]) if ii != 1 else sznl[1] for ii in range(len(sznl))]
         vmax = max([np.abs(sz).max().values for sz in [sznl[0], sznl[2]]])
      else:
         vmax = max([np.abs(sz).max().values for sz in sznl])
         vmin = min([sz.min().values for sz in sznl])
         zero_centered = True
         if vmax - vmin < 0.25 * vmax:
            zero_centered = False

      plt.rc('font', size=16)
      plt.rcParams['figure.figsize'] = (30, 12)
      contourfkwargs = {'cmap': 'coolwarm' if zero_centered else 'rainbow', 'norm': colors.CenteredNorm(vcenter=0, halfrange=vmax) if zero_centered else colors.Normalize(vmin=vmin, vmax=vmax), 'levels': 11}
      subplot_kw = dict(xlim=(-1, 1), ylim=(0, 1))
      fig, axes = plt.subplots(2, 5, layout='constrained', sharey=True, subplot_kw=subplot_kw)
      plt.suptitle(f"{dv} [{dss[0][dv].attrs['units']}]\n{dss[0][dv].attrs['long_name']}\n{dss[0][dv].attrs['cell_methods']}")

      for ii, sz in enumerate(['JJA', 'SON']):
         for jj, cs in enumerate(CASES):
            ax, csf = axes[ii, jj], None
            if DO_DIFF:
               if jj == 1:
                  csf = ax.contourf(YSCL(sznl[jj][ymeth[0]]), ZSCL(sznl[jj][zmeth[0]]), sznl[jj].sel(season=sz).data, levels=11, cmap='rainbow')
               else:
                  csf = ax.contourf(YSCL(sznl[jj][ymeth[0]]), ZSCL(sznl[jj][zmeth[0]]), sznl[jj].sel(season=sz).data, levels=contourfkwargs['levels'], cmap='bwr', norm=colors.TwoSlopeNorm(0))
            else:
               csf = ax.contourf(YSCL(sznl[jj][ymeth[0]]), ZSCL(sznl[jj][zmeth[0]]), sznl[jj].sel(season=sz).data, **contourfkwargs)
            if ii == 0 and jj == 0:
               ax.set_ylabel('Depth [m]')
               ax.set_yticks(ZLOC, ZLAB)
               ax.invert_yaxis()
               if not DO_DIFF:
                  contourfkwargs['levels'] = csf.levels * 1.25
            plt.colorbar(csf, ax=ax)
            ax.set_xticks(YLOC, YLAB)
            ax.set_xlabel('Latitude [°]')

      plt.savefig(os.path.join(DIRO, '%s.png' % dv), bbox_inches='tight')
      plt.close()

   outds.to_netcdf(os.path.join(DIRO, 'sznlzm.nc'))

def main_plot():
   ds = xr.open_dataset(os.path.join(DIRO, 'sznlzm.nc'))
   print(ds)

   for dv in ds.data_vars:
      print('Plotting', str(dv))
      plt.rc('font', size=16)
      plt.rcParams['figure.figsize'] = (30, 12)
      #contourfkwargs = {'cmap': 'coolwarm' if zero_centered else 'rainbow', 'norm': colors.CenteredNorm(vcenter=0, halfrange=vmax) if zero_centered else colors.Normalize(vmin=vmin, vmax=vmax), 'levels': 11}
      subplot_kw = dict(xlim=(-1, 1), ylim=(0, 1))
      fig, axes = plt.subplots(2, 3, layout='constrained', sharey=True, subplot_kw=subplot_kw)
      #plt.suptitle(f"{dv} [{dss[0][dv].attrs['units']}]\n{dss[0][dv].attrs['long_name']}\n{dss[0][dv].attrs['cell_methods']}")
      plt.suptitle(str(dv))

      #print(ds[dv], ds[dv].dims[-1])
      ydim = ds[dv].dims[-1]
      levels = (11, 11)
      if str(dv) in CLEVS:
         levels = (np.arange(*CLEVS[str(dv)][0]), np.arange(*CLEVS[str(dv)][1]))
      for ii, sz in enumerate(ds['season']):
         for jj, cs in enumerate(ds['case']):
            ax, csf = axes[ii, jj], None
            ctrlnorm = colors.TwoSlopeNorm(0) if str(dv) in ZEROCTR else None
            if DO_DIFF:
               if jj == 1:
                  csf = ax.contourf(YSCL(ds[ydim]), ZSCL(ds['zl']), ds[dv].sel(case=cs, season=sz).data, levels=levels[0], cmap='PiYG', norm=ctrlnorm)
               else:
                  ctrl = ds[dv].isel(case=1).sel(season=sz).data
                  csf = ax.contourf(YSCL(ds[ydim]), ZSCL(ds['zl']), ds[dv].sel(case=cs, season=sz).data - ctrl, levels=levels[1], cmap='bwr', norm=colors.TwoSlopeNorm(0), extend='both')
            else:
               csf = ax.contourf(YSCL(ds[ydim]), ZSCL(ds['zl']), ds[dv].sel(case=cs, season=sz).data, levels=levels[0], cmap='PiYG', norm=ctrlnorm)
            if ii == 0 and jj == 0:
               ax.set_ylabel('Depth [m]')
               ax.set_yticks(ZLOC, ZLAB)
               ax.invert_yaxis()

               if str(dv) in LIMDEPTH:
                  ax.set_yticks(LIMLOC, LIMLAB)
                  ax.set_ylim(ZSCL(800), ZSCL(0)) 

            plt.colorbar(csf, ax=ax)
            ax.set_xticks(YLOC, YLAB)
            ax.set_xlabel('Latitude [°]')
   
      plt.savefig(os.path.join(DIRO, '%s.png' % str(dv)), bbox_inches='tight')
      plt.close()

if __name__ == '__main__':
   if mode == 'compute':
      main()
   if mode == 'plot':
      main_plot()
   if mode not in MODES:
      raise ValueError ('Invalid mode for script')
