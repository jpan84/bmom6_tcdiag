import numpy as np
import xarray as xr
from sznl_funcs import stack_hemi_sznl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gspec

TCDENS = '/glade/u/home/jpan/aquaptc/tempest/250725_density_sznl/tcdens.nc'
OCN_YZ = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/ocn_yz_plts_diff/sznlzm.nc'
H1IALL = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/tcfields2mps_%s_full/means_all.nc'
H1ITCS = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/tcfields2mps_%s_full/means_tcs.nc'

ALIS = ['250415_unseed', '250417_ctrl', '250416_seed1x1']

ZSCL = np.vectorize(lambda zl: zl / 2850 if zl <= 2000 else 2000 / 2850 + (zl - 2000) / 2000 * 850 / 2850)
ZLAB = np.arange(1000, 5000, 1000)
ZLOC = ZSCL(ZLAB)

LIMLAB = np.arange(200, 1500, 200)
LIMLOC = ZSCL(LIMLAB)

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)

def main():
   tcdens = xr.open_dataset(TCDENS)
   ocn_yz = xr.open_dataset(OCN_YZ)
   h1iall = xr.concat([xr.open_dataset(H1IALL % al).expand_dims(case=[al]) for al in ALIS], dim='case')
   h1itcs = xr.concat([xr.open_dataset(H1ITCS % al).expand_dims(case=[al]) for al in ALIS], dim='case')

   #mirror and stack the seasons of TC density
   for dv in tcdens.data_vars:
      stacked = stack_hemi_sznl(tcdens[dv], antisym=False, latnm='lat')
      tcdens = tcdens.assign(variables={dv: stacked})
   tcdens = tcdens.sel(season=['JJA', 'SON'])

   print(tcdens)
   print(ocn_yz)
   print(h1iall)
   print(h1itcs)

   #TAUX and ACE density above vmo
   sf_cfkwargs = {'cmap': 'coolwarm', 'levels': 2.**np.arange(-2, 7, 1), 'norm': colors.SymLogNorm(2.**-1)}
   sf_cfkwargs['levels'] = np.concatenate((-sf_cfkwargs['levels'][::-1], sf_cfkwargs['levels']))
   sf_ctkwargs = {'colors': 'black', 'levels': 2.**np.arange(-2, 7, 1)}
   sf_ctkwargs['levels'] = np.concatenate((-sf_ctkwargs['levels'][::-1], sf_ctkwargs['levels']))
   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (24, 9)
   fig, axes = plt.subplots(2, 3, layout='constrained')
   for ii, szn in enumerate(ocn_yz['season']):
      for jj in range(len(ALIS)):
         outer_spec = axes[ii, jj].get_subplotspec()
         gs = gspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_spec, height_ratios=[1, 3])
         ax_top = fig.add_subplot(gs[0])
         ax_bot = fig.add_subplot(gs[1])

         sf = (ocn_yz['vhGM'] + ocn_yz['vmo']).sel(season=szn).isel(case=jj)
         ax_bot.contourf(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / 1e10, **sf_cfkwargs)
         ax_bot.contour(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / 1e10, **sf_ctkwargs)

   fig.tight_layout()
   plt.show()


if __name__ == '__main__':
   main()
