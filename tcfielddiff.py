import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

CTLDIR = './tcfields4mps_250417_ctrl/'
EXPDIR = './tcfields4mps_250416_seed1x1/'
EXPDIR = './tcfields4mps_250415_unseed/'
#EXPDIR = './tcfields4mps_250417_ctrl/'
DIRO = './tcfieldsdiff/'

#FLD = ['hflso_neg', 'hflso_pos']
#SGN = -1
latvar = 'yh'
#FLD = ['PRECT']
#SGN = 1
latvar = 'latitudes'
FLD = ['TAUX_neg', 'TAUX_pos']
SGN = 1
latvar = 'latitudes'
FLD = ['VQ850_neg', 'VQ850_pos']
SGN = 1
latvar = 'latitudes'

lv = 2.501e6 #default from MOM6
rho_l = 1.0e3 #default from CAM zm micro
mps2mmph = 1e3 * 3.6e3
mps2mmpd = 1e3 * 8.64e4

AVGLATBND = (20, 30)

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-50, 51, 10)
MT = YSCL(np.array([10, 20]))
STZ = YSCL(np.array([20, 30]))

wgtavg = lambda var, wgt: (var * wgt).sum() / wgt.sum()

################### namelist-like section #########################################
#TAUX UNSEED
LTR = '(a)'
TTL = 'TAUX UNSEED$-$CTRL'
ZLIMS = (-.02, .02)
ZTICKS = None
NONTC = False
PLTFRAC = False
ISCTRL = False
LGND = False
MULTC = 1
OUTNM = 'UNSEED-CTRL_TAUX.png'

#TAUX CTRL
LTR = '(b)'
TTL = 'TAUX CTRL [N m$^{-2}$]'
ZLIMS = (-.12, .12)
ZTICKS = np.arange(-.12, .13, .02)
NONTC = False
PLTFRAC = True
ISCTRL = True
LGND = True
MULTC = 5
OUTNM = 'CTRL_TAUX.png'

#TAUX SEED
LTR = '(c)'
TTL = 'TAUX SEED$-$CTRL'
ZLIMS = (-.02, .02)
ZTICKS = None
NONTC = False
PLTFRAC = False
ISCTRL = False
LGND = False
MULTC = 1
OUTNM = 'SEED-CTRL_TAUX.png'

#VQ850 UNSEED
LTR = '(a)'
TTL = 'UNSEED$-$CTRL $v$\'$q$\'$_{850}$'
ZLIMS = (-3.5e-4, 3.5e-4)
ZTICKS = None
NONTC = False
PLTFRAC = False
ISCTRL = False
LGND = False
MULTC = 1
OUTNM = 'UNSEED-CTRL_VQ850.png'

#VQ850 CTRL
#LTR = '(b)'
#TTL = 'CTRL $v$\'$q$\'$_{850}$ [m s$^{-1}$ kg kg$^{-1}$]'
#ZLIMS = (-7e-3, 7e-3)
#ZTICKS = None
#NONTC = False
#PLTFRAC = True
#ISCTRL = True
#LGND = True
#MULTC = 5
#OUTNM = 'CTRL_VQ850.png'

#VQ850 SEED
#LTR = '(c)'
#TTL = 'SEED$-$CTRL $v$\'$q$\'$_{850}$'
#ZLIMS = (-7e-3, 7e-3)
#ZTICKS = None
#NONTC = False
#PLTFRAC = False
#ISCTRL = False
#LGND = False
#MULTC = 1
#OUTNM = 'SEED-CTRL_VQ850.png'


###################################################################################

if not os.path.exists(DIRO):
   os.makedirs(DIRO)

ctlall = xr.open_dataset(os.path.join(CTLDIR, 'means_all.nc'))
ctltcs = xr.open_dataset(os.path.join(CTLDIR, 'means_tcs.nc'))
expall = xr.open_dataset(os.path.join(EXPDIR, 'means_all.nc'))
exptcs = xr.open_dataset(os.path.join(EXPDIR, 'means_tcs.nc'))
dss = [ctlall, expall, ctltcs, exptcs]

#P minus E
#dss = [ds.interp(coords=dict(latitudes=ds.yh.data)) for ds in dss]

flds = [SGN * sum([ds[dv] for dv in FLD]) for ds in dss]
#diff = [flds[1] - flds[0], flds[-1] - flds[-2]]

#P minus E
#flds = [((ds['hflso_neg'] + ds['hflso_pos']) / lv / rho_l + ds['PRECT'].rename(latitudes='yh')) * mps2mmpd for ds in dss]
ctlfac = 0 if ISCTRL else 1
diff = [flds[1] - ctlfac*flds[0], flds[-1] - ctlfac*flds[-2]]
#flds = [-((ds['hflso_neg'] + ds['hflso_pos']) / lv / rho_l) * mps2mmpd for ds in dss] #E
#flds = [ds['PRECT'] * mps2mmpd for ds in dss] #P
#diff = [flds[1] - flds[0], flds[-1] - flds[-2], (flds[1] - flds[-1]) - (flds[0] - flds[-2])]

# !outside Portion of field outside TCs instead of all
#diff = [(flds[1] - flds[-1]) - (flds[0] - flds[-2]), flds[-1] - flds[-2]]

plt.rc('font', size=14)
sinlat = YSCL(ctlall[latvar])
ax0 = plt.axes()
ax1 = plt.gca().twinx() if PLTFRAC else None
if PLTFRAC:
   ax0.set_ylabel('TC fractional contribution')
   if FLD[0] == 'VQ850_neg':
      ax1.plot(sinlat.sel(latitudes=slice(14, 40)), (diff[1].sum('season') / diff[0].sum('season')).sel(latitudes=slice(14, 40)), linewidth=3, color='black') #ratio for ctrl only
      ax1.set_ylim(-0.1, 0.1)
   elif FLD[0] == 'TAUX_neg':
      ax1.plot(sinlat.sel(latitudes=slice(18, 33)), (diff[1].sum('season') / diff[0].sum('season')).sel(latitudes=slice(18, 33)), linewidth=3, color='black') #ratio for ctrl only
      ax1.set_ylim(-.06, .06)
   else:
      ax1.plot(sinlat, diff[1].sum('season') / diff[0].sum('season'), linewidth=2.5, color='black') #ratio for ctrl only
   
for ii, szn in enumerate(diff[0]['season']):
   ax0.plot(sinlat, diff[0].sel(season=szn), color=['blue', 'orange'][ii], label=str(szn.data) + ' total')
   ax0.plot(sinlat, MULTC * diff[1].sel(season=szn), color=['blue', 'orange'][ii], linestyle='dashed', label=str(szn.data) + ' TCs')
   if NONTC:
      ax0.plot(sinlat, diff[2].sel(season=szn), color=['blue', 'orange'][ii], linestyle='dotted', label=str(szn.data) + ' non-TC') # !outside
   ax0.set_xlim(YSCL(YLAB[0]), YSCL(YLAB[-1]))
   ax0.set_xticks(YSCL(YLAB), YLAB)
   ax0.axhline(0, linewidth=0.5, color='gray')
   ax0.axvspan(*MT, fc='purple', alpha=.08)
   ax0.axvspan(*STZ, fc='yellow', alpha=.12)
   [plt.axvline(ll, linewidth=0.5, color='gray') for ll in YSCL(YLAB)]

   toavg_all = diff[0].sel(indexers={'season': szn, latvar: slice(*AVGLATBND)})
   toavg_tcs = diff[1].sel(indexers={'season': szn, latvar: slice(*AVGLATBND)})
   coslat = np.cos(np.deg2rad(toavg_all[latvar]))

   print(str(szn.data), 'average over', AVGLATBND)
   print('All:', wgtavg(toavg_all, coslat).data)
   print('TCs:', wgtavg(toavg_tcs, coslat).data)

   #plt.ylim(-18, 32)
   #plt.yticks(np.arange(-16, 33, 4))
   #plt.title('LHFLX SEED–CTRL')
   #plt.ylim(-.02, .02)
   #plt.title('TAUX UNSEED–CTRL')

   #plt.ylim(-1.3 / 5, 1.3 / 5)
   #plt.title('P [mm d$^{-1}$] CTRL')
   plt.title(LTR, loc='left')
   plt.title(TTL)
   ax0.set_ylim(*ZLIMS)
   if ZTICKS is not None:
      ax0.set_yticks(ZTICKS)
   #ax0.set_ylim(0, 18)
   #ax1.set_ylim(0, 0.18)

   #ax0.set_ylim(-.12, .12)
   #ax0.set_yticks(np.arange(-.12, .13, .02))
   #plt.title('TAUX CTRL [N m$^{-2}$]')

   # !outside
   #plt.ylim(-35, 35)
   #plt.yticks(np.arange(-32, 33, 4))
   #plt.title('LHFLX SEED–CTRL')

   #plt.ylim(-2.1, 2.1)
   #plt.title('P–E [mm d$^{-1}$] CTRL')

   #ax0.set_ylim(-7e-3, 7e-3)
   #ax0.set_ylim(-3e-4, 3e-4)
   #ax1.set_ylim(-.12, .12)
   #plt.title('UNSEED$-$CTRL $v$\'$q$\'$_{850}$ [m s$^{-1}$ kg kg$^{-1}$]')

   ax0.set_xlabel('Latitude')
   if LGND:
      legend = ax0.legend(ncol=2, fontsize=11, framealpha=1)
      legend.set_zorder(9)
   #ax1.legend(*ax0.get_legend_handles_labels(), fontsize=11, framealpha=1, ncol=3)

plt.savefig(os.path.join(DIRO, OUTNM))
plt.show()
