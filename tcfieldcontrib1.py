#compute the fraction or component of an h1i field that is accounted for by TCs
#Use raw h1i output and NodeFileFilter'ed h1i output, and yhour eddy fields

import sys
import os
import math
import operator
import uxarray as ux
import xarray as xr
import numpy as np
from sznl_funcs import monthly2sznl, stack_hemi_sznl
import pltsettings
import matplotlib.pyplot as plt

MODES = ['compute', 'plot']
mode = sys.argv[1]

ARCHV = '/glade/campaign/univ/upsu0032/jpan_tcfields/'
ALIAS = '250417_ctrl'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s' % ALIAS
DIRO = './tcfields2mps_%s/' % ALIAS
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

RAWALL = 'hist_0010_h1i/cat_h1i.nc'
RAWTCS = 'TC_R2_masked/*.h1i.*.nc'
EDDALL = 'yhoureddy/yhoureddy_h1i.nc'
EDDTCS = 'yhoureddy_TC_R2_masked/yhoureddy_h1i.nc'

g = 9.81
unsigned_vars = ['PRECT']
signed_vars = dict(TAUX=(-1, 'TAUX', False), TAUY=(-1, 'TAUY', True), UMF500=(-1/g, 'OMEGA500', False), UMF850=(-1/g, 'OMEGA850', False))#template: (...take the product of these scalars/data_vars..., antisym?)

eddy_fluxes=dict(VTBOT=('VBOT', 'TBOT', True), VT850=('V850', 'T850', True), VT500=('V500', 'T500', True), VT200=('V200', 'T200', True),\
                 VUBOT=('VBOT', 'UBOT', True), VU850=('V850', 'U850', True), VU500=('V500', 'U500', True), VU200=('V200', 'U200', True),\
                 VQ850=('V850', 'Q850', True), VZ850=('V850', 'Z850', True), VZ500=('V500', 'Z500', True),\
                 WT850=(-1, 'OMEGA850', 'T850', False), WT500=(-1, 'OMEGA500', 'T500', False),\
                 WU850=(-1, 'OMEGA850', 'U850', False), WU500=(-1, 'OMEGA500', 'U500', False),\
                 WQ850=(-1, 'OMEGA850', 'Q850', False), WZ850=(-1, 'OMEGA850', 'Z850', False), WZ500=(-1, 'OMEGA500', 'Z500', False)) 
eddy_fluxes=dict(VT850=('V850', 'T850', True),\
                 VU200=('V200', 'U200', True),\
                 VQ850=('V850', 'Q850', True),\
                 WT850=(-1, 'OMEGA850', 'T850', False),\
                 WU500=(-1, 'OMEGA500', 'U500', False),\
                 WQ850=(-1, 'OMEGA850', 'Q850', False))

zmlats = (-90, 90.1, 0.5)
LATLAB = np.arange(-90, 91, 30)

ZMLATS = np.arange(*zmlats)
SINLAT = np.sin(np.deg2rad(ZMLATS))
LSTYS = ['solid', 'dashed']
LCLRS = ['blue', 'orange']
NEWLATS = np.arange(-50, 51, 10)

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)
   pltsettings.set1()

   print('Opening datasets...')
   rawallds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, RAWALL))
   rawtcsds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, RAWTCS))
   eddallds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, EDDALL))
   eddtcsds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, EDDTCS))

   outallds, outtcsds = compute_unsigned_flds(rawallds, rawtcsds, unsigned_vars)

   sgnallds, sgntcsds = compute_signed_flds(rawallds, rawtcsds, signed_vars)
   outallds, outtcsds = xr.merge([outallds, sgnallds]), xr.merge([outtcsds, sgntcsds])

   eddallds, eddtcsds = compute_signed_flds(eddallds, eddtcsds, eddy_fluxes)
   outallds, outtcsds = xr.merge([outallds, eddallds]), xr.merge([outtcsds, eddtcsds])

   print('Saving .nc outputs...')
   outallds.to_netcdf(os.path.join(DIRO, 'means_all.nc'))
   outtcsds.to_netcdf(os.path.join(DIRO, 'means_tcs.nc'))

   print(sys.argv[0], 'done.')

def main_plot():
   allds = xr.open_dataset(os.path.join(DIRO, 'means_all.nc'))
   tcsds = xr.open_dataset(os.path.join(DIRO, 'means_tcs.nc'))

   plt.rc('font', size=16)

   for dv in allds.data_vars:
      print('Plotting', dv)
      plt.rcParams['figure.figsize'] = (10, 12)
      subplot_kw = dict(xlim=(-0.75, 0.75), sharex=True)
      fig, axes = plt.subplots(2, 1, subplot_kw=subplot_kw, gridspec_kw=dict(height_ratios=[3, 1]))
      axes[0].hlines(0, -1, 1, colors='black', linestyles='dotted')
      ratio = tcsds[dv] / allds[dv]
      for tt, szn in enumerate(allds['season']):
         axes[1].plot(SINLAT, ratio[tt], color=LCLRS[tt])
         axes[1].set_title('TCs / all')
         for ii, pltda in enumerate([allds[dv], tcsds[dv]]):
            axes[0].plot(SINLAT, pltda.sel(season=szn), label=str(szn.values), color=LCLRS[tt], linestyle=LSTYS[ii])
            axes[0].set_title(str(dv))
      axes[0].legend()
      [ax.set_xticks(np.sin(np.deg2rad(NEWLATS)), NEWLATS) for ax in axes]
      axes[1].set_xlabel('Latitude [°]')
      fig.tight_layout()
      plt.savefig(os.path.join(DIRO, str(dv)))
      plt.close()

#compute the product of terms (scalar or var name) in template
def template_prod(ds, templ):
   terms = [ds[var] if type(var) == str else float(var) for var in templ]
   return math.prod(terms)

#Compute seasonal and zonal mean after all physical quantity computation/masking is complete
def all_and_TC_to_sznlzm(allda, tcsda, antisym=False):
   monzm = [da.groupby('time.month').mean().zonal_mean(lat=ZMLATS) for da in [allda, tcsda]]
   return [stack_hemi_sznl(monthly2sznl(da), antisym=antisym) for da in monzm]

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
            retall = xr.Dataset(data_vars={fl + sgn[0]: sznzm[0].to_xarray()})
            rettcs = xr.Dataset(data_vars={fl + sgn[0]: sznzm[1].to_xarray()})
         else:
            retall = retall.assign(variables={fl + sgn[0]: sznzm[0].to_xarray()})
            rettcs = rettcs.assign(variables={fl + sgn[0]: sznzm[1].to_xarray()})
   return retall, rettcs


if __name__ == '__main__':
   if mode == 'compute':
      main()
   if mode == 'plot':
      main_plot()
   if mode not in MODES:
      raise ValueError('invalid mode for script')
