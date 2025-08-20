#compute the fraction or component of an h1i field that is accounted for by TCs
#Use raw h1i output and NodeFileFilter'ed h1i output, and yhour eddy fields

import sys
import os
import math
import operator
import uxarray as ux
import xarray as xr
from dask.diagnostics import ProgressBar
import numpy as np
from sznl_funcs import monthly2sznl, stack_hemi_sznl
import pltsettings
import matplotlib.pyplot as plt

MODES = ['compute', 'plot']
mode = sys.argv[1]

SPTH = 4

#ARCHV = '/glade/campaign/univ/upsu0032/jpan_tcfields/'
ARCHV = '/glade/derecho/scratch/jpan/jpan_tcfields/'
ALIAS = '250415_unseed'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s' % ALIAS
DIRO = './tcfields%dmps_%s/' % (SPTH, ALIAS)
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

RAWALL = 'hist_0012-0014_h1i/cat_h1i_0012-0014.nc'
RAWTCS = 'hist_0012-0014_h1i/TC_R%d_masked_h1i_0012-0014.nc' % SPTH
EDDALL = 'hist_0012-0014_h1i/yhoureddy_h1i_0012-0014.nc'
EDDTCS = 'hist_0012-0014_h1i/TC_R%d_masked_yhoureddy_h1i_0012-0014.nc' % SPTH

OCNALL = 'hist_0012-0014_h1i/cat_sfc_0012-0014.nc'
OCNTCS = 'hist_0012-0014_h1i/TC_R%d_masked_sfc_0012-0014.nc' % SPTH
OCAALL = 'hist_0012-0014_h1i/yhoureddy_sfc_0012-0014.nc'
OCATCS = 'hist_0012-0014_h1i/TC_R%d_masked_yhoureddy_sfc_0012-0014.nc' % SPTH

g = 9.81

#CAM fields
unsigned_vars = ['PRECT']
signed_vars = dict(TAUX=(-1, 'TAUX', False), TAUY=(-1, 'TAUY', True), UMF500=(-1/g, 'OMEGA500', False), UMF850=(-1/g, 'OMEGA850', False))#template: (...take the product of these scalars/data_vars..., antisym?)

eddy_fluxes=dict(VTBOT=('VBOT', 'TBOT', True), VT850=('V850', 'T850', True), VT500=('V500', 'T500', True), VT200=('V200', 'T200', True),\
                 VUBOT=('VBOT', 'UBOT', True), VU850=('V850', 'U850', True), VU500=('V500', 'U500', True), VU200=('V200', 'U200', True),\
                 VQ850=('V850', 'Q850', True), VZ850=('V850', 'Z850', True), VZ500=('V500', 'Z500', True),\
                 WT850=(-1, 'OMEGA850', 'T850', False), WT500=(-1, 'OMEGA500', 'T500', False),\
                 WU850=(-1, 'OMEGA850', 'U850', False), WU500=(-1, 'OMEGA500', 'U500', False),\
                 WQ850=(-1, 'OMEGA850', 'Q850', False), WZ850=(-1, 'OMEGA850', 'Z850', False), WZ500=(-1, 'OMEGA500', 'Z500', False)) 
###eddy_fluxes=dict(VT850=('V850', 'T850', True),\
###                 VU200=('V200', 'U200', True),\
###                 VQ850=('V850', 'Q850', True),\
###                 WT850=(-1, 'OMEGA850', 'T850', False),\
###                 WU500=(-1, 'OMEGA500', 'U500', False),\
###                 WQ850=(-1, 'OMEGA850', 'Q850', False))
###eddy_fluxes, signed_vars = dict(VT850=('V850', 'T850', True)), dict(TAUY=(-1, 'TAUY', True))

#MOM fields
ocn_signed = dict(hflso=(1, 'hflso', False), hfsso=(1, 'hfsso', False))
ocn_anom = dict(rlntds=(1, 'rlntds', False), rsntds=(1, 'rsntds', False), hflso_a=(1, 'hflso', False), hfsso_a=(1, 'hfsso', False))

zmlats = (-90, 90.1, 0.5)
LATLAB = np.arange(-90, 91, 30)

ZMLATS = np.arange(*zmlats)
SINLAT = np.sin(np.deg2rad(ZMLATS))
LSTYS = ['solid', 'dashed']
LCLRS = ['blue', 'orange']
NEWLATS = np.arange(-50, 51, 10)

RATLIMS = dict(PRECT=0.7, TAUX_neg=0.5, TAUX_pos=1, TAUY_neg=0.6, TAUY_pos=1, UMF500_neg=0.35, UMF500_pos=0.65,\
            VQ850_neg=0.6, VQ850_pos=0.6, VT850_neg=0.6, VT850_pos=0.6, VU200_neg=0.3, VU200_pos=0.25, WQ850_neg=0.5, WQ850_pos=0.5, WT850_neg=0.4, WT850_pos=0.6, WU500_neg=0.8, WU500_pos=0.7)
YLIMS = dict(PRECT=(0, 2.2e-7), TAUX_neg=(-0.12, .01), TAUX_pos=(-5e-3, .03), TAUY_neg=(-.1, .1), TAUY_pos=(-.1, .1), UMF500_neg=(-.012, .012), UMF500_pos=(-.012, .012),\
            VQ850_neg=(-.009, .009), VQ850_pos=(-.009, .009), VT850_neg=(-5, 5), VT850_pos=(-5, 5), VU200_neg=(-100, 100), VU200_pos=(-100, 100), WQ850_neg=(-1e-4, 1e-5), WQ850_pos=(-2e-5, 5e-4),\
            WT850_neg=(-0.15, 0), WT850_pos=(0, .15), WU500_neg=(-.8, .05), WU500_pos=(-.05, .8))

RLIMS_NET = dict(TAUX_net=0.55, TAUY_net=0.65, UMF500_net=0.5, VQ850_net=0.5, VT850_net=0.6, VU200_net=0.25, VU850_net=0.65, VZ500_net=0.65, WQ850_net=0.5, WT850_net=0.5, WU500_net=0.7, WZ850_net=0.75)
YLIMS_NET = dict(TAUX_net=(-.13, .2), TAUY_net=(-.05, .05), UMF500_net=(-3e-3, 8e-3), VQ850_net=(-8e-3, 8e-3), VT850_net=(-18, 15), VU200_net=(-85, 85), VU850_net=(-15, 12), VZ500_net=(-20, 20),\
                WQ850_net=(-1e-5, 4e-4), WT850_net=(-2e-2, .18), WU500_net=(-0.12, 0.35), WZ850_net=(-3.6, 0.1))

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)
   pltsettings.set1()

   print('Opening datasets...')
   rawallds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, RAWALL))
   rawtcsds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, RAWTCS))
   eddallds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, EDDALL))
   eddtcsds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, EDDTCS))
   ocnallds = xr.open_dataset(os.path.join(ARCHV, CASE, OCNALL))
   ocntcsds = xr.open_dataset(os.path.join(ARCHV, CASE, OCNTCS))
   ocaallds = xr.open_dataset(os.path.join(ARCHV, CASE, OCAALL))
   ocatcsds = xr.open_dataset(os.path.join(ARCHV, CASE, OCATCS))

   print('Setting up CAM fields...')
   outallds, outtcsds = compute_unsigned_flds(rawallds, rawtcsds, unsigned_vars)

   sgnallds, sgntcsds = compute_signed_flds(rawallds, rawtcsds, signed_vars)
   outallds, outtcsds = xr.merge([outallds, sgnallds]), xr.merge([outtcsds, sgntcsds])

   eddallds, eddtcsds = compute_signed_flds(eddallds, eddtcsds, eddy_fluxes)
   outallds, outtcsds = xr.merge([outallds, eddallds]), xr.merge([outtcsds, eddtcsds])

   print('Setting up MOM fields...')
   ocnallds, ocntcsds = compute_signed_flds(ocnallds, ocntcsds, ocn_signed)
   ocaallds, ocatcsds = compute_signed_flds(ocaallds, ocatcsds, ocn_anom)
   outallds, outtcsds = xr.merge([outallds, ocnallds, ocaallds]), xr.merge([outtcsds, ocntcsds, ocatcsds])

   print('Saving .nc outputs...')
   with ProgressBar():
      outallds.to_netcdf(os.path.join(DIRO, 'means_all.nc'))
   with ProgressBar():
      outtcsds.to_netcdf(os.path.join(DIRO, 'means_tcs.nc'))

   print(sys.argv[0], 'done.')

def main_plot():
   allds = xr.open_dataset(os.path.join(DIRO, 'means_all.nc'))
   tcsds = xr.open_dataset(os.path.join(DIRO, 'means_tcs.nc'))

   plt.rc('font', size=16)

   for dv in allds.data_vars:
      print('Plotting', dv)
      #print(allds[dv].dims)
      plt.rcParams['figure.figsize'] = (10, 12)
      subplot_kw = dict(xlim=(-0.75, 0.75), sharex=True)
      fig, axes = plt.subplots(2, 1, subplot_kw=subplot_kw, gridspec_kw=dict(height_ratios=[3, 1]))
      axes[0].hlines(0, -1, 1, colors='black', linestyles='dotted')
      ratio = (tcsds[dv] / allds[dv]).clip(min=0)
      sinlats = SINLAT
      for tt, szn in enumerate(allds['season']):
         if 'yh' in ratio.dims:
            sinlats = np.sin(np.deg2rad(ratio['yh']))
         axes[1].plot(sinlats, ratio[tt], color=LCLRS[tt])
         axes[1].set_title('TCs / all')
         for ii, pltda in enumerate([allds[dv], tcsds[dv]]):
            axes[0].plot(sinlats, pltda.sel(season=szn), label=str(szn.values), color=LCLRS[tt], linestyle=LSTYS[ii])
            axes[0].set_title(str(dv))
      axes[0].legend()
      [ax.set_xticks(np.sin(np.deg2rad(NEWLATS)), NEWLATS) for ax in axes]
      axes[1].set_xlabel('Latitude [°]')
      if str(dv) in RATLIMS:
         axes[1].set_ylim(0, RATLIMS[dv])
      if str(dv) in YLIMS:
         axes[0].set_ylim(*YLIMS[dv])
      fig.tight_layout()
      plt.savefig(os.path.join(DIRO, str(dv)))
      plt.close()

      if '_pos' in str(dv):
         totname = str(dv).replace('_pos', '_net')
         negname = str(dv).replace('_pos', '_neg')
         #print(str(dv), totname)
         plt.rcParams['figure.figsize'] = (10, 12)
         subplot_kw = dict(xlim=(-0.75, 0.75), sharex=True)
         fig, axes = plt.subplots(2, 1, subplot_kw=subplot_kw, gridspec_kw=dict(height_ratios=[3, 1]))
         axes[0].hlines(0, -1, 1, colors='black', linestyles='dotted')
         netall, nettcs = allds[dv] + allds[negname], tcsds[dv] + tcsds[negname]
         #ratio = nettcs / netall #raw ratio is unstable at zero crossings
         posrat, negrat = tcsds[dv] / allds[dv], tcsds[negname] / allds[negname]
         ratio = (np.abs(allds[dv]) * posrat + np.abs(allds[negname]) * negrat)\
                    / (np.abs(allds[dv]) + np.abs(allds[negname]))
         ratio = ratio.clip(min=0)
         for tt, szn in enumerate(allds['season']):
            axes[1].plot(sinlats, ratio[tt], color=LCLRS[tt])
            axes[1].set_title('TCs / all')
            for ii, pltda in enumerate([netall, nettcs]):
               axes[0].plot(sinlats, pltda.sel(season=szn), label=str(szn.values), color=LCLRS[tt], linestyle=LSTYS[ii])
               axes[0].set_title(totname)
         axes[0].legend()
         [ax.set_xticks(np.sin(np.deg2rad(NEWLATS)), NEWLATS) for ax in axes]
         axes[1].set_xlabel('Latitude [°]')
         if totname in RLIMS_NET:
            axes[1].set_ylim(0, RLIMS_NET[totname])
         if totname in YLIMS_NET:
            axes[0].set_ylim(*YLIMS_NET[totname])
         fig.tight_layout()
         plt.savefig(os.path.join(DIRO, totname))
         plt.close()

#compute the product of terms (scalar or var name) in template
def template_prod(ds, templ):
   terms = [ds[var] if type(var) == str else float(var) for var in templ]
   return math.prod(terms)

#Compute seasonal and zonal mean after all physical quantity computation/masking is complete
def all_and_TC_to_sznlzm(allda, tcsda, antisym=False):
   typ = type(allda)
   if typ == ux.UxDataArray:
      monzm = [da.groupby('time.month').mean().zonal_mean(lat=ZMLATS) for da in [allda, tcsda]]
      return [stack_hemi_sznl(monthly2sznl(da), antisym=antisym) for da in monzm]
   if typ == xr.DataArray: #hard-coded for MOM6 tracer points
      monzm = [da.groupby('time.month').mean().mean(dim='xh') for da in [allda, tcsda]]
      return [stack_hemi_sznl(monthly2sznl(da), antisym=antisym, latnm='yh') for da in monzm]

def compute_unsigned_flds(allds, tcsds, flds):
   retall, rettcs = None, None
   for fl in flds:
      print('Processing %s unsigned...' % fl)
      sznzm = all_and_TC_to_sznlzm(allds[fl], tcsds[fl])
      if retall is None:
         retall = xr.Dataset(data_vars={fl: sznzm[0].to_xarray()})
         rettcs = xr.Dataset(data_vars={fl: sznzm[1].to_xarray()})
      else:
         retall = retall.assign(variables={fl: sznzm[0].to_xarray()})
         rettcs = rettcs.assign(variables={fl: sznzm[1].to_xarray()})
   return retall, rettcs

def compute_signed_flds(allds, tcsds, flds_dict):
   retall, rettcs = None, None
   for fl in flds_dict:
      print('Processing %s signed...' % fl)
      antisym = flds_dict[fl][-1]
      templ = flds_dict[fl][:-1]
      allfld = template_prod(allds, templ)
      tcsfld = template_prod(tcsds, templ)
      sgnlat = np.sign(allds['lat']) if antisym else 1

      for sgn in [('_pos', operator.gt), ('_neg', operator.lt)]:
         allsgn = (sgn[1](allfld * sgnlat, 0)) * allfld
         tcssgn = (sgn[1](tcsfld * sgnlat, 0)) * tcsfld
         sznzm = all_and_TC_to_sznlzm(allsgn, tcssgn, antisym=antisym)
         if retall is None:
            retall = xr.Dataset(data_vars={fl + sgn[0]: sznzm[0]})#.to_xarray()})
            rettcs = xr.Dataset(data_vars={fl + sgn[0]: sznzm[1]})#.to_xarray()})
         else:
            retall = retall.assign(variables={fl + sgn[0]: sznzm[0]})#.to_xarray()})
            rettcs = rettcs.assign(variables={fl + sgn[0]: sznzm[1]})#.to_xarray()})
   return retall, rettcs


if __name__ == '__main__':
   if mode == 'compute':
      main()
   if mode == 'plot':
      main_plot()
   if mode not in MODES:
      raise ValueError('invalid mode for script')
