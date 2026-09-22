import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post')
from paths import ARCHRT, ALIA, CTLIX, CASENAMES, CAMGR, IXHORS
import consts as c
import uxarray as ux
import xarray as xr
import numpy as np
from sznl_funcs import stack_hemi_sznl, monthly2sznl, agg_time

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker

h0 = '/glade/campaign/univ/upsu0032/jpan_aquaptc/%s/atm/hist/*.h0a.001[0-9]*.nc'
FILO = 'MMC_native_0010-0019.nc'

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-90, 91, 10)
YLOC = YSCL(YLAB)

def main_compute():
   dss = [ux.open_mfdataset(CAMGR, h0 % cs).expand_dims(case=[ALIA[ii]]) for ii, cs in enumerate(CASENAMES)]
   for ii in range(1, len(dss)):
      dss[ii].uxgrid = dss[0].uxgrid
   ds = ux.concat(dss, dim='case')

   print('Setting up coords...')
   aterm = ds['hyai'] * c.P0
   bterm = ds['hybi'] * ds['PS']
   p_ilev = aterm + bterm
   dp3d = p_ilev.diff('ilev').rename(dict(ilev='lev')).assign_coords(lev=ds['lev'])
   ds = ds.assign(variables=dict(dp3d=dp3d, p_ilev=p_ilev, coslat=np.cos(np.deg2rad(ds['lat']))))

   vdp = ds['V'] * ds['dp3d']
   vdp_zm = vdp.zonal_mean((-90, 90, 1.5)).assign_coords(case=ds['case'], time=ds['time'])
   p_ilev_zm = ds['p_ilev'].zonal_mean((-90, 90, 1.5)) #the integrals values correspond to interfaces

   #cum = vdp_zm.cumulative('lev').sum()
   vdp_zm_ilev = vdp_zm.rename(lev='ilev')
   
   vdp_zm_ilev = xr.concat(
       [xr.zeros_like(vdp_zm_ilev.isel(ilev=0)), vdp_zm_ilev],
       dim='ilev',
   )
   
   vdp_zm_ilev = vdp_zm_ilev.assign_coords(ilev=p_ilev_zm.ilev)
   
   cum = vdp_zm_ilev.cumsum('ilev')
   p_tgt = np.concatenate((ds['lev'].data, [1010.])) * 100.

   myint = lambda data, pi: np.interp(p_tgt, pi, data, left=np.nan, right=np.nan)

   cum = vdp_zm_ilev.cumsum('ilev').chunk({'ilev': -1})
   p_ilev_zm = p_ilev_zm.chunk({'ilev': -1})
   #print(p_ilev_zm.isel(ilev=-1).max().values, p_ilev_zm.isel(ilev=-1).min().values)

   cum_p = xr.apply_ufunc(
           myint, cum, p_ilev_zm,
           input_core_dims=[['ilev'], ['ilev']],
           output_core_dims=[['plev']],
           vectorize=True,
           dask='parallelized',
           output_dtypes=[cum.dtype],
           dask_gufunc_kwargs={'output_sizes': {'plev': len(p_tgt)}},
   ).assign_coords(plev=p_tgt)

   print(cum_p)

   mmc_sf = cum_p * 2 * np.pi * c.a_e * np.cos(np.deg2rad(cum_p['latitudes'])) / c.g
   mmc_sf.to_dataset(name='MMC_SF').to_netcdf(FILO)


def main_plot():
   ds = xr.open_dataset(FILO)

   #test_plot = ds['MMC_SF'].mean(dim=['time', 'case'])

   #plt.contourf(test_plot['latitudes'], test_plot['plev'], test_plot.T)
   #plt.colorbar()
   #plt.show()

   #hy = ds['MMC_SF'].mean(dim='time') #test placeholder
   hy = agg_time(ds['MMC_SF'], antisym=True)

   print('Plotting streamfunctions...')
   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (30, 12)
   contourfkwargs = {'cmap': 'bwr', 'levels': 2.**np.arange(-2, 7, 1), 'norm': colors.SymLogNorm(2.**-1)} #coolwarm for diff
   contourfkwargs['levels'] = np.concatenate((-contourfkwargs['levels'][::-1], contourfkwargs['levels']))
   contourkwargs = {'colors': 'black', 'levels': 2.**np.arange(-2, 7, 1)}
   contourkwargs['levels'] = np.concatenate((-contourkwargs['levels'][::-1], contourkwargs['levels']))
   clabelkwargs = {'inline': 1, 'fontsize': 10, 'colors': 'black', 'fmt': '%.1f'}
   subplot_kw = dict(xlim=(-1, 1), ylim=(100, 1000), yscale='log')
   fig, axes = plt.subplots(2, 3, layout='constrained', sharey=True, subplot_kw=subplot_kw)

   hy /= 1e10

   for ii, ax in enumerate(axes.ravel()):
      ix = IXHORS[ii]
      toplt = hy.isel(case=ix)
      contourfkwargs['cmap'] = 'coolwarm'
      if ix != CTLIX:
         toplt = (hy.isel(case=ix) - hy.isel(case=CTLIX)) * 10
         contourfkwargs['cmap'] = 'bwr'

      CSF = ax.contourf(np.sin(np.deg2rad(hy['latitudes'])), hy['plev'] / 100., toplt.data.T, **contourfkwargs)
      CS1 = ax.contour(np.sin(np.deg2rad(hy['latitudes'])), hy['plev'] / 100., hy.isel(case=ix).data.T, **contourkwargs) 
      if ii == 0:
         ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
         ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
         ax.set_ylabel('Pressure [hPa]')
         ax.set_yticks(np.arange(100, 1001, 100))
         cb = plt.colorbar(CSF, ax=axes)
         cbt = cb.get_ticks()
         cbt = np.concatenate((cbt[:cbt.size // 2 + 1], -cbt[:cbt.size // 2 + 1][::-1]))
         cb.set_ticks(cbt)
         cb.set_ticklabels([(('%d' % t) if abs(t) >= 1 else t) for t in cbt])
         ax.invert_yaxis()
      ax.set_xticks(YLOC, labels=[yl if abs(yl) <= 60 else '' for yl in YLAB])

      #ax.set_title('%s (10$^{%d}$ kg s$^{-1}$)' % (plottitles[sfi], expo))
      ax.set_xlabel('Latitude [°]')

   #fig.tight_layout()
   plt.savefig(FILO + '.svg', bbox_inches='tight')
   plt.show()


if __name__ == '__main__':
   if len(sys.argv) == 1:
      main_compute()
   elif sys.argv[1] == 'plot':
      main_plot()
