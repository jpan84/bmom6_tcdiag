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
TAUX = [(1, 'TAUX'), ('mean', 'lon')]
DIAB = [(cp / g, 'DTCOND'), (cp / g, 'QRL'), (cp / g, 'QRS'), ('mean', 'lon')]

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
   ###if DO_DIFF:
   ###   oht[0], oht[2] = oht[0] - oht[1], oht[2] - oht[1]
   siny_q = np.sin(np.deg2rad(oht[0]['yq']))
   ###print(sinlat_q.values) #TODO: always remove the last (northernmost) yq point to make hemi sym

   print('\tPlotting OHT...')
   plt_cases(siny_q, oht, 'OHT', '[W]', linestyle='solid')

   print('\nComputing AHT...')
   msefl = [derive_flx(ds, MSEFLX, 'lat', True) for ds in camdss]
   ####if DO_DIFF:
   #print(camdss[0]['lat'])
   lat_h = msefl[0]['lat']
   #print(lat_h)
   latcirc = 2 * np.pi * a * np.cos(np.deg2rad(lat_h))
   aht = [latcirc / g * fl.integrate('plev') for fl in msefl]
   sinlat_h = np.sin(np.deg2rad(lat_h))

   print('\tPlotting AHT...')
   plt_cases(sinlat_h, aht, 'AHT', '[W]')

   print('\nSumming AHT + OHT...')
   aht_q = [da.interp(coords=dict(lat=oht[0]['yq'].data)) for da in aht]
   aoht = [oht[ii] + ada.data.T for ii, ada in enumerate(aht_q)]
   print('\tPlotting AOHT...')
   plt_cases(sinlat_q, aoht, 'AOHT', '[W]')

   print('\nComputing O Angular Momentum Transport...')
   vo_h = [ds['vo'].rename(dict(yq='yh')).interp(coords=dict(yh=ds['yh'])) for ds in momdss] #non-conservative interp
   momdss = [ds.assign(variables=dict(vo_h=vo_h[ii])) for ii, ds in enumerate(momdss)]
   #print(momdss[0]['vo_h'])
   #print(momdss[0]['rhoinsitu'])
   r_ax_yh = a * np.cos(np.deg2rad(momdss[0]['yh']))
   #print(r_ax)
   #print(momdss[0]['uv'])
   oamf_r = [r_ax_yh**2 * derive_flx(ds, OamFLX_rot, 'yh', True) for ds in momdss]
   oamf_u = [r_ax_yh * derive_flx(ds, OamFLX_u, 'yh', True) for ds in momdss]
   oamt = [2 * np.pi * r_ax_yh * (oamf_u[ii] + oamf_r[ii]).integrate('zl') for ii in range(len(oamf_r))]
   print('\tPlotting OAMT...')
   plt_cases(sinlat_h, oamt, 'OAMT', '[N m]')

   print('\nComputing A Angular Momentum Transport...')
   aamf_r = [r_ax**2 * derive_flx(ds, AamFLX_rot, 'lat', True) for ds in camdss]
   aamf_u = [r_ax * derive_flx(ds, AamFLX_u, 'lat', True) for ds in camdss]
   aamt = [latcirc / g * (aamf_r[ii] + aamf_u[ii]).integrate('plev') for ii in range(len(aamf_r))]
   print('\tPlotting AAMT...')
   plt_cases(sinlat_h, aamt, 'AAMT', '[N m]')

   print('\nSumming AAMT + OAMT...')
   AOamT = [oamt[ii] + ada.sel(lat=oamt[ii]['yh'].data).data.T for ii, ada in enumerate(aamt)]
   print('\tPlotting AOamT...')
   plt_cases(sinlat_h, AOamT, 'AO_AMT', '[N m]')

   print('\nComputing AHU...')
   ahu = [derive_flx(ds, AHU, 'lat', False) for ds in camdss]
   print('\tPlotting AHU...')
   plt_cases(sinlat_h, ahu, 'AHU', '[W m$^{-2}$]')

   print('\nComputing sfcHU...')
   shu = [derive_flx(ds, OHU_CAM, 'lat', False) for ds in camdss]
   print('\tPlotting sfcHU...')
   plt_cases(sinlat_h, shu, 'sfcHU', '[W m$^{-2}$]')

   print('\nComputing RESTOM...')
   rom = [derive_flx(ds, RESTOM, 'lat', False) for ds in camdss]
   print('\tPlotting RESTOM...')
   plt_cases(sinlat_h, rom, 'RESTOM', '[W m$^{-2}$]')

   print('\nComputing TAUX torque...')
   taux_tor = [r_ax * derive_flx(ds, TAUX, 'lat', False) for ds in camdss]
   print('\tPlotting TAUX torque...')
   plt_cases(sinlat_h, taux_tor, 'atmo_sfc_torque', '[N m m$^{-2}$]')

   print('\nComputing ocean stress torques...')
   otor = [r_ax * derive_flx(ds, OMU, 'yh', False) for ds in momdss]
   print('\tPlotting ocean stress torques...')
   plt_cases(sinlat_h, otor, 'ocn_stress_torque', '[N m m$^{-2}$]')

   print('\nComputing diabatic heating...')
   diab = [derive_flx(ds, DIAB, 'lat', False).integrate('plev') for ds in camdss]
   print('\tPlotting diabatic heating...')
   plt_sznl(sinlat_h, diab, 'diabatic_heating', '[W m$^{-2}$]')
   

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
   plt.rcParams['figure.figsize'] = (30, 6)
   subplot_kw = dict(xlim=(-1, 1), sharey=(not DO_DIFF))
   fig, axes = plt.subplots(1, 3, layout='constrained', subplot_kw=subplot_kw)
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
