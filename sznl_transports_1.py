import os
import math
import numpy as np
import xarray as xr
#from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt

from sznl_funcs import monthly2sznl, stack_hemi_sznl
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

RESTOM = [(1, 'FSNT'), (-1, 'FLNT'), ('mean', 'lon')]
AHU = [(1, 'FSNT'), (-1, 'FSNS'), (1, 'FLNS'), (-1, 'FLNT'), (1, 'SHFLX'), (1, 'LHFLX'), ('mean', 'lon')]
OHU_CAM = [(1, 'FSNS'), (-1, 'FLNS'), (-1, 'SHFLX'), (-1, 'LHFLX'), ('mean', 'lon')]
MSEFLX = [(g, 'V', 'Z3'), (cp, 'VT'), (lv, 'VQ'), ('mean', 'lon')]


OHU_MOM = [(1, 'rsntds'), (1, 'rlntds'), (1, 'hfsso'), (1, 'hflso'), ('mean', 'xh')]
OHT = [(1, 'T_ady_2d'), (1, 'T_diffy_2d'), ('sum', 'xh')]
OMU = [(1, 'taux_bot'), (1, 'tauuo'), ('mean', 'xq')]

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   print('Opening history files...')
   camdss = [xr.open_mfdataset(CAMH % cs) for cs in CASES]
   momdss = [xr.open_mfdataset(MOMH % cs) for cs in CASES]

   print('\nComputing OHT...')
   oht = [derive_flx(ds, OHT) for ds in momdss]
   sinlat_q = np.sin(np.deg2rad(oht[0]['yq']))

   print('\tPlotting...')
   plt.rcParams['figure.figsize'] = (30, 6)
   subplot_kw = dict(xlim=(-1, 1), sharey=(not DO_DIFF))
   fig, axes = plt.subplots(1, 3, subplot_kw=subplot_kw)
   for ii, pltda in enumerate(oht):
      plt_sznl(sinlat_q, pltda, 'OHT', '[W]', title=ALIASES[ii], linestyle='solid', ax=aces[ii])
   plt.savefig(os.path.join(DIRO, 'OHT.png'), bbox_inches='tight')
   plt.close()

#Take the product of terms in template tuples, then sum over tuples
#Final element in template gives zonal aggregation method
def derive_flx(ds, template):
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

   return sznl

###   exit()
###   ################################################
###
###   print('Opening history files...')
###   momdss = [xr.open_mfdataset(os.path.join(ARCHV, CASE1, MOMH))]
###   camdss = [xr.open_mfdataset(os.path.join(ARCHV, CASE1, CAMH))]
###   if CASE2 is None:
###      momdss.append(xr.zeros_like(momdss[0]))
###      camdss.append(xr.zeros_like(camdss[0]))
###   else:
###      momdss.append(xr.open_mfdataset(os.path.join(ARCHV, CASE2, MOMH)))
###      camdss.append(xr.open_mfdataset(os.path.join(ARCHV, CASE2, CAMH)))
###   camdss = [ds.rename(dict(lat='lat_h')) for ds in camdss]
###   momdss = [ds.rename(dict(yh='lat_h')) for ds in momdss]
###
###   print('Computing OHT...')
###   oht = [monmeans((ds.T_ady_2d + ds.T_diffy_2d).sum(dim='xh')) for ds in momdss]
###   oht = oht[0] - oht[1]
###   oht_sznl = wmat @ oht.data
###   #plt.plot(oht.yq, oht.mean(dim='time').values)
###   #plt.show()
###   sinlat_q = np.sin(np.deg2rad(oht.yq))
###   plt_sznl(sinlat_q, oht_sznl, 'OHT', '[W]')
###
###   print('Computing AHT...')
###   gpf = [g * monzm(ds['V']) * monzm(ds['Z3']) for ds in camdss]
###   shf = [cp * monzm(ds['VT']) for ds in camdss]
###   lhf = [lv * monzm(ds['VQ']) for ds in camdss]
###   mseflx = sum([f[0] - f[1] for f in [gpf, shf, lhf]])
###   latcirc = 2 * np.pi * a * np.cos(np.deg2rad(mseflx.lat_h))
###   aht = latcirc / g * mseflx.integrate('plev')
###   aht = aht.transpose('month', ...)
###   aht_q = aht.interp(coords=dict(lat_h=oht.yq.data)).rename(dict(lat_h='yq'))
###   print(aht.shape)
###   aht_sznl = wmat @ aht.data
###   aht_sznl_q = wmat @ aht_q.data
###   sinlat_cam_h = np.sin(np.deg2rad(mseflx.lat_h))
###   plt_sznl(sinlat_cam_h, aht_sznl, 'AHT', '[W]')
###
###   print('Plotting AHT + OHT...')
###   plt_sznl(sinlat_q, oht_sznl + aht_sznl_q, 'AOHT', '[W]')
###
###   print('Computing AMT...')
###   umf = [monzm(ds['VU']) for ds in camdss]
###   umf = umf[0] - umf[1]
###   amt = latcirc / g * umf.integrate('plev')
###   amt = amt.transpose('month', ...)
###   amt_sznl = wmat @ amt.data
###   plt_sznl(sinlat_cam_h, amt_sznl, 'AMT', '[N]')
###
###   print('Computing OMT...')
###   omf = [monmeans(ds['rhoinsitu'] * ds['uv']).mean(dim='xh') for ds in momdss]
###   omf = omf[0] - omf[1]
###   omt = latcirc * omf.integrate('zl')
###   omt = omt.transpose('month', ...)
###   omt_sznl = wmat @ omt.data
###   sinlat_h = np.sin(np.deg2rad(omt.lat_h))
###   plt_sznl(sinlat_h, omt_sznl, 'OMT', '[N]')
###
###   print('Computing OMU...')
###   omu = [monmeans((ds['tauuo'] + ds['taux_bot']).mean(dim='xq')) for ds in momdss]
###   omu = omu[0] - omu[1]
###   omu_sznl = wmat @ omu.data
###   plt_sznl(sinlat_q, omu_sznl, 'OMU', '[N m$^{-2}$]')
###
###   print('Computing TAUX...')
###   taux = [monzm(ds['TAUX']) for ds in camdss]
###   taux = taux[0] - taux[1]
###   taux_sznl = wmat @ taux.transpose('month', ...).data
###   plt_sznl(sinlat_cam_h, taux_sznl, 'TAUX', '[N m$^{-2}$]')
###
###   print('Computing AHU...')
###   ahu = [sum([monzm(tm[0] * ds[tm[1]]) for tm in formulae.AHU]) for ds in camdss]
###   ahu = ahu[0] - ahu[1]
###   ahu = ahu.transpose('month', ...)
###   ahu_sznl = wmat @ ahu.data
###   plt_sznl(sinlat_cam_h, ahu_sznl, 'AHU', '[W m$^{-2}$]')
###
###   print('Computing sfc heat uptake...')
###   shu = [sum([monzm(tm[0] * ds[tm[1]]) for tm in formulae.OHU_CAM]) for ds in camdss]
###   shu = shu[0] - shu[1]
###   shu = shu.transpose('month', ...)
###   shu_sznl = wmat @ shu.data
###   plt_sznl(sinlat_cam_h, shu_sznl, 'sfcHU', '[W m$^{-2}$]')
###
###   print('Computing RESTOM...')
###   rom = [sum([monzm(tm[0] * ds[tm[1]]) for tm in formulae.RESTOM]) for ds in camdss]
###   rom = rom[0] - rom[1]
###   rom = rom.transpose('month', ...)
###   rom_sznl = wmat @ rom.data
###   plt_sznl(sinlat_cam_h, rom_sznl, 'RESTOM', '[W m$^{-2}$]')
###
###   print('Computing diabatic heating...')
###   ttnd_diab = [monzm(ds['DTCOND'] + ds['QRL'] + ds['QRS']) for ds in camdss]
###   diab_int = cp / g * (ttnd_diab[0] - ttnd_diab[1]).integrate('plev')
###   diab_int = wmat @ diab_int.transpose('month', ...).data
###   plt_sznl(sinlat_cam_h, diab_int, 'integ_Qdot', '[W m$^{-2}$]')
###
###   print('Done.')
###
####vals should have shape (time, lat)
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
###
###   if ax is None:
###      ax = plt.axes()
###   for tt in range(vals.shape[0]):
###      plt.plot(sinlats, vals[tt, :], label=SZNS[tt], color=lncolors[tt], linestyle=linestyle)
###   if True or not CASE2 is None:
###      plt.hlines(0, -1, 1, color='black', linestyle='dotted')
###   plt.xlabel('Lat [°]')
###   plt.ylabel(name + ' ' + units)
###   if units == '[W]' and CASE2 is None:
###      plt.ylim(-2e16, 2e16)
###   plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=10))
###   plt.xlim(-1, 1)
###   plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
###   plt.savefig(os.path.join(DIRO, '%s.%s.png' % (CASE1.split('.')[-1], name)), bbox_inches='tight')
###   if close:
###      plt.close()

if __name__ == '__main__':
   main()
