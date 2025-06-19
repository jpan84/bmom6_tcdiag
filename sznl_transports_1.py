import os
import math
import numpy as np
import xarray as xr
#from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt

import sznl_funcs
import rad_sfc_aht_oht as formulae

DIRO = 'sznl_transports_new'
ARCHV = '/glade/derecho/scratch/jpan/archive/'
CAMH = '/glade/derecho/scratch/jpan/archive/%s/atm/hist_regrid_mom_onpres/*.h0a.*.nc'
MOMH = '/glade/derecho/scratch/jpan/archive/%s/ocn/hist/*mom6.hm*[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'
CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
ALIASES = ['UNSEED', 'CTRL', 'SEED']
DO_DIFF = False

LATLAB = np.array([-90., -60., -30., 0., 30., 60., 90.])
lncolors = ['blue', 'orange']
monmeans = lambda da: da.groupby('time.month').mean()

g = 9.81
cp = 1004
lv = 2500840
a = 6.371e6
OM = 7.29e-5

RESTOM = [(1, 'FSNT'), (-1, 'FLNT'), ('mean', 'lon')]
AHU = [(1, 'FSNT'), (-1, 'FSNS'), (1, 'FLNS'), (-1, 'FLNT'), (1, 'SHFLX'), (1, 'LHFLX'), ('mean', 'lon')]
OHU_CAM = [(1, 'FSNS'), (-1, 'FLNS'), (-1, 'SHFLX'), (-1, 'LHFLX'), ('mean', 'lon')]
MSEFLX = [(g, 'V', 'Z3'), (cp, 'VT'), (lv, 'VQ'), ('mean', 'lon')]


OHU_MOM = [(1, 'rsntds'), (1, 'rlntds'), (1, 'hfsso'), (1, 'hflso'), ('mean', 'xh')]
OHT = [(1, 'T_ady_2d'), (1, 'T_diffy_2d'), ('sum', 'xh')]
OMU = [(1, 'taux_bot'), (1, 'tauuo'), ('mean', 'xq')]

OamFLX_rot = [(OM, 'rhoinsitu', 'vo_h'), ('mean', 'xh')] #vo must be interpolated from yq:point to yh:mean
OamFLX_u = [(1, 'rhoinsitu', 'uv'), ('mean', 'xh')]
AamFLX_rot = [(OM, 'V'), ('mean', 'lon')]
AamFLX_u = [(1, 'VU'), ('mean', 'lon')]

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   print('Opening history files...')
   camdss = [xr.open_mfdataset(CAMH % cs) for cs in CASES]
   momdss = [xr.open_mfdataset(MOMH % cs) for cs in CASES]

   print('\nComputing OHT...')
   oht = [derive_flx(ds, OHT, 'yq', True).isel(yq=slice(0, -1)) for ds in momdss]
   if DO_DIFF:
      oht[0], oht[2] = oht[0] - oht[1], oht[2] - oht[1]
   sinlat_q = np.sin(np.deg2rad(oht[0]['yq']))
   #print(sinlat_q.values) #TODO: always remove the last (northernmost) yq point to make hemi sym

   print('\tPlotting OHT...')
   plt_cases(sinlat_q, oht, 'OHT', '[W]', linestyle='solid')

   print('\nComputing AHT...')
   msefl = [derive_flx(ds, MSEFLX, 'lat', True) for ds in camdss]
   #if DO_DIFF:
   lat_h = msefl[0]['lat']
   latcirc = 2 * np.pi * a * np.cos(np.deg2rad(lat_h))
   aht = [latcirc / g * fl.integrate('plev') for fl in msefl]
   sinlat_h = np.sin(np.deg2rad(lat_h))

   print('\tPlotting AHT...')
   plt_cases(sinlat_h, aht, 'AHT', '[W]')

   print('\nSumming AHT + OHT...')
   aht_q = [aht.interp(coords=dict(lat=oht[0]['yq'].data))]
   aoht = [oht[ii] + ada.data for ii, ada in enumerate(aht_q)]
   print('\tPlotting AOHT...')
   plt_cases(sinlat_q, aoht, 'AOHT', '[W]')

   print('\nComputing O Angular Momentum Transport...')
   vo_h = [ds['vo'].rename(dict(yq='yh')).interp(coords=dict(yh=ds['yh'])) for ds in momdss] #non-conservative interp
   momdss = [ds.assign(variables=dict(vo_h=vo_h[ii])) for ii, ds in enumerate(momdss)]
   r_ax = latcirc / 2 / np.pi
   oamf_r = [r_ax**2 * derive_flx(ds, OamFLX_rot, 'yh', True) for ds in momdss]
   oamf_u = [r_ax * derive_flx(ds, OamFLX_u, 'yh', True) for ds in momdss]
   oamt = [latcirc * (oamf_u[ii] + oamf_r[ii]).integrate('zl') for ii in range(len(oamf_r))]
   print('\tPlotting OAMT...')
   plt_cases(sinlat_h, oamt, 'OAMT', '[N m]')

   print('\nComputing A Angular Momentum Transport...')
   aamf_r = [r_ax**2 * derive_flx(ds, AamFLX_rot, 'lat', True) for ds in camdss]
   aamf_u = [r_ax * derive_flx(ds, AamFLX_u, 'lat', True) for ds in camdss]
   aamt = [latcirc / g * (aamf_r[ii] + aamf_u[ii]).integrate('plev') for ii in range(len(aamf_r))]
   print('\tPlotting AAMT...')
   plt_cases(sinlat_h, aamt, 'AAMT', '[N m]')

   print('\nSumming AAMT + OAMT...')
   AOamT = [oamt[ii] + ada.data for ii, ada in enumerate(aamt)]
   print('\tPlotting AOamT...')
   plt_cases(sinlat_h, AOamT, 'AO_AMT', '[N m]')

#Take the product of terms in template tuples, then sum over tuples
#Final element in template gives zonal aggregation method
def derive_flx(ds, template, latnm, antisym):
   termtups = template[:-1]
   lonagg = template[-1]

   rawvars = [[ds[var] if type(var) == str else float(var) for var in tup] for tup in termtups]
   prods = [math.prod(rv) for rv in rawvars]

   mthly = monmeans(sum(prods))
   monzm = None
   if lonagg[0] == 'mean':
      monzm = mthly.mean(dim=lonagg[1])
   elif lonagg[0] == 'sum':
      monzm = mthly.sum(dim=lonagg[1])
   sznl = sznl_funcs.monthly2sznl(monzm)
   stacked = sznl_funcs.stack_hemi_sznl(sznl, antisym=antisym, latnm=latnm) #TODO: check yh and yq symmetry

   return stacked


def plt_cases(sinlat, das, *args, **kwargs):
   plt.rcParams['figure.figsize'] = (20, 4)
   subplot_kw = dict(xlim=(-1, 1), sharey=(not DO_DIFF))
   fig, axes = plt.subplots(1, 3, subplot_kw=subplot_kw)
   for ii, pltda in enumerate(das):
      plt_sznl(sinlat, pltda, *args, **kwargs, title=ALIASES[ii], ax=axes[ii])
   plt.savefig(os.path.join(DIRO, '%s.png' % args[0]), bbox_inches='tight')
   plt.close()

###pltda should have shape (season, lat)
def plt_sznl(sinlats, pltda, name, units, title='', linestyle='solid', ax=None):
   if ax is None:
      ax = plt.axes()
   for tt, szn in enumerate(pltda['season']):
      ax.plot(sinlats, pltda.sel(season=szn), label=str(szn.values), color=lncolors[tt])
   ax.set_xlabel('Lat [°]')
   ax.set_ylabel(name + ' ' + units)
   ax.set_title(title)
   ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=12))
   ax.set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
   ax.hlines(0, -1, 1, color='black', linestyle='dotted')

if __name__ == '__main__':
   main()
