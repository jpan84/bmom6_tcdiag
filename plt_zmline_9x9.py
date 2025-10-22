import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os

DIRI = './linevslat_h0a_diff'
FILI = 'sznlzm.nc'
TTLS = ['UNSEED$-$CTRL', 'CTRL', 'SEED$-$CTRL']
PLTVARS = ['SHFLX', 'LHFLX', 'FLUS']
lncolors = ['blue', 'orange']

YLIMS = dict(SHFLX=[(-.5, .5), (-18, 18), (-2.7, 2.7)], LHFLX=[(-3.2, 3.2), (-10, 210), (-22, 22)], FLUS=[(-0.8, 1.2), (240, 520), (-8, 12)])
VARLABS = dict(SHFLX='SH [W m$^{-2}$]', LHFLX='LH [W m$^{-2}$]', FLUS='$F_{LW}^↑$ [W m$^{-2}$]')
YSCL = lambda deglat: np.sin(np.deg2rad(deglat))
YLAB = np.arange(-60, 61, 10)
YLOC = YSCL(YLAB)

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
      ax.set_xticks(YLOC, ['' if yl % 30 else yl for yl in YLAB])
      [ax.axvline(yl, c='gray', lw=0.5) for yl in YLOC]
      ax.axvspan(-1, YLOC[0], fill=True, fc='gray', zorder=2, alpha=0.6)
      ax.axvspan(YLOC[-1], 1, fill=True, fc='gray', zorder=2, alpha=0.6)

      ax.set_ylim(*YLIMS[dv][jj])
      ax.set_ylabel(VARLABS[dv], fontsize=16)
      ax.set_title(mosaic_labels[ii][jj], loc='left')
      ax.legend(ncol=2, fontsize=14)

fig.tight_layout()
plt.savefig(os.path.join(DIRI, '_'.join(PLTVARS) + '.svg'))
#plt.show()
exit()



print('\tPlotting...')
plt.rcParams['figure.figsize'] = (24, 18)
subplot_kw = dict(xlim=(-1, 1), sharey=(not DO_DIFF))
fig, axes = plt.subplots(1, 3, subplot_kw=subplot_kw)

for ii, pltda in enumerate(sznzm):
   ax = axes[ii]
   for tt, szn in enumerate(pltda['season']):
      ax.plot(sinlat, pltda.sel(season=szn), label=str(szn.values), color=lncolors[tt])
   #plt.plot(sinlat, line1.values)
   #if str(dv) == 'TS' and not DO_DIFF:
   #   plt.hlines(273.15 + 26.5, -1, 1, colors='red', linestyles='dashed')
   if DO_DIFF and (ii == 0 or ii == 2):
      ax.hlines(0, -1, 1, colors='black', linestyles='dashed')
   ax.set_xlabel('Lat [°]')
   try:
      ax.set_ylabel(str(dv) + ' [' + str(das[1].units) + ']')
   except:
      ax.set_ylabel(dv)
   ax.set_title(ALIASES[ii] + (' difference' if (DO_DIFF and (ii == 0 or ii == 2)) else '') )
   ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=12))
   ax.set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
   [ax.axvline(np.sin(np.deg2rad(ll)), c='gray', lw=0.5) for ll in LATLAB]

if DO_DIFF:
   ###make the difference panels share y
   #ylims = [ax.get_ylim() for ax in axes]
   #ylims.pop(1)
   #ylims = np.array(ylims)
   #miny, maxy = ylims.min(), ylims.max()
   #axes[0].set_ylim(miny, maxy)
   #axes[2].set_ylim(miny, maxy)
   for aa in [axes[0], axes[2]]:
      ylims = aa.get_ylim()
      maxy = max(np.abs(ylims))
      aa.set_ylim(-maxy, maxy)
plt.savefig(os.path.join(OUTDIR, '%s_sznl.png' % dv), bbox_inches='tight')
plt.close()
