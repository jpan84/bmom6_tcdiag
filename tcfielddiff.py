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
#latvar = 'yh'
FLD = ['TAUX_neg', 'TAUX_pos']
SGN = 1
latvar = 'latitudes'
#FLD = ['VQ850_neg', 'VQ850_pos']
#SGN = 1
#latvar = 'latitudes'

lv = 2.501e6 #default from MOM6
rho_l = 1.0e3 #default from CAM zm micro
mps2mmph = 1e3 * 3.6e3
mps2mmpd = 1e3 * 8.64e4

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-50, 51, 10)

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
diff = [flds[1] - flds[0], flds[-1] - flds[-2]]

#P minus E
#flds = [((ds['hflso_neg'] + ds['hflso_pos']) / lv / rho_l + ds['PRECT'].rename(latitudes='yh')) * mps2mmpd for ds in dss]
#diff = [flds[1] - 0*flds[0], flds[-1] - 0*flds[-2]]

# !outside Portion of field outside TCs instead of all
#diff = [(flds[1] - flds[-1]) - (flds[0] - flds[-2]), flds[-1] - flds[-2]]

plt.rc('font', size=14)
sinlat = YSCL(ctlall[latvar])
for ii, szn in enumerate(diff[0]['season']):
   plt.plot(sinlat, diff[0].sel(season=szn), color=['blue', 'orange'][ii], label=str(szn.data) + ' total')
   plt.plot(sinlat, diff[1].sel(season=szn), color=['blue', 'orange'][ii], linestyle='dashed', label=str(szn.data) + ' TCs')
   #plt.plot(sinlat, diff[0].sel(season=szn), color=['blue', 'orange'][ii], linestyle='dotted', label=str(szn.data) + '_non-TC') # !outside
   plt.xlim(YSCL(YLAB[0]), YSCL(YLAB[-1]))
   plt.xticks(YSCL(YLAB), YLAB)
   plt.axhline(0, linewidth=0.5, color='gray')
   [plt.axvline(ll, linewidth=0.5, color='gray') for ll in YSCL(YLAB)]

   #plt.ylim(-18, 32)
   #plt.yticks(np.arange(-16, 33, 4))
   #plt.title('LHFLX UNSEED–CTRL')
   plt.ylim(-.02, .02)
   plt.title('TAUX UNSEED–CTRL')

   #plt.ylim(-.12, .12)
   #plt.title('TAUX CTRL')

   # !outside
   #plt.ylim(-35, 35)
   #plt.yticks(np.arange(-32, 33, 4))
   #plt.title('LHFLX SEED–CTRL')

   #plt.ylim(-2.1, 2.1)
   #plt.title('P–E [mm d$^{-1}$] CTRL')

   #plt.ylim(-7e-3, 7e-3)
   #plt.title('V\'Q\'850 [m/s kg/kg] UNSEED$-$CTRL')

   #plt.legend()

plt.savefig(os.path.join(DIRO, 'UNSEED-CTRL_TAUX.png'))
plt.show()
