import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os

DIRI = './linevslat_h0a_diff'
FILI = 'sznlzm.nc'
TTLS = ['UNSEED$-$CTRL', 'CTRL', 'SEED$-$CTRL']
PLTVARS = ['SHFLX', 'LHFLX', 'FLUS']
#PLTVARS = ['FSNS', 'FLNS', 'FLNSC']
lncolors = ['blue', 'orange']

YLIMS = dict(SHFLX=[(-.5, .5), (-18, 18), (-2.7, 2.7)], LHFLX=[(-3.2, 3.2), (-10, 210), (-22, 22)], FLUS=[(-0.8, 1.2), (240, 520), (-8, 12)],\
              FSNS=[(-1.6, 1.2), (0, 280), (-24, 18)], FLNS=[(-1.6, 1.2), (0, 102), (-24, 18)], FLNSC=[(-1.6, 1.2), (0, 102), (-24, 18)])
VARLABS = dict(SHFLX='SH [W m$^{-2}$]', LHFLX='LH [W m$^{-2}$]', FLUS='$F_{LW}^↑$ [W m$^{-2}$]', FSNS='$F_{SW}^{net}$ [W m$^{-2}$]',\
                FLNS='$F_{LW}^{net}$ [W m$^{-2}$]', FLNSC='$F_{LW}^{net,clear}$ [W m$^{-2}$]')
YSCL = lambda deglat: np.sin(np.deg2rad(deglat))
YLAB = np.arange(-60, 61, 10)
YLOC = YSCL(YLAB)
MT = YSCL(np.array([10, 20]))
STZ = YSCL(np.array([20, 30]))

ds = xr.open_dataset(os.path.join(DIRI, FILI))
sinlat = YSCL(ds['latitudes'])

plt.rc('font', size=18)
mosaic_labels_flat = [f'({chr(i)})' for i in range(97, 97 + 9)]
mosaic_labels = [mosaic_labels_flat[i:i + 3] for i in range(0, 9, 3)]
fig, axes = plt.subplots(len(mosaic_labels), len(mosaic_labels[0]), figsize=(16, 12), sharex=True, sharey=False, subplot_kw=dict(xlim=(-1, 1)))

for ii, dv in enumerate(PLTVARS):
   for jj, cs in enumerate(ds['case']):
      ax = axes[ii][jj]
      for tt, szn in enumerate(ds['season']):
         ax.plot(sinlat, ds[dv].sel(case=cs, season=szn), label=str(szn.values), color=lncolors[tt])

      if ii == 0:
         ax.set_title(TTLS[jj])

      ax.axhline(0, c='gray', lw=0.5)
      ax.tick_params(right=True, labelbottom=True)
      ax.set_xlabel('Latitude [°]', fontsize=15)
      ax.set_xticks(YLOC, ['' if yl % 20 else yl for yl in YLAB])
      [ax.axvline(yl, c='gray', lw=0.5) for yl in YLOC]
      ax.axvspan(-1, YLOC[0], fill=True, fc='gray', zorder=2, alpha=0.77)
      ax.axvspan(YLOC[-1], 1, fill=True, fc='gray', zorder=2, alpha=0.77)
      ax.axvspan(*MT, fill=True, fc='purple', zorder=2, alpha=0.08)
      ax.axvspan(*STZ, fill=True, fc='yellow', zorder=2, alpha=0.12)

      ax.set_ylim(*YLIMS[dv][jj])
      ax.set_ylabel(VARLABS[dv], fontsize=16)
      ax.set_title(mosaic_labels[ii][jj], loc='left')
      ax.legend(ncol=2, fontsize=14)

fig.tight_layout()
plt.savefig(os.path.join(DIRI, '_'.join(PLTVARS) + '.svg'))
plt.show()
