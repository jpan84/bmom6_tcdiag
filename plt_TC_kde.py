import pickle
import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cftime

DIRI = '/glade/u/home/jpan/aquaptc/tempest/'
FILI = ['250415_unseed.kde.pickle', '250417_ctrl.kde.pickle', '250416_seed1x1.kde.pickle']
PQTS = ['250415_unseed.parquet', '250417_ctrl.parquet', '250416_seed1x1.parquet']
SSTS = 'zm_sst_ydaymean.nc'
TTLS = [('(a)', 'UNSEED$-$CTRL'), ('(b)', 'CTRL'), ('(c)', 'SEED$-$CTRL')]

TICKDATES = [cftime.datetime(1, ii, 1, calendar='noleap') for ii in range(1, 13)]
TICKDOY = [int(dt.strftime('%-j')) for dt in TICKDATES]

pkls = []
for fl in FILI:
   with open(os.path.join(DIRI, fl), 'rb') as ofl:
      pkls.append(pickle.load(ofl))
#print(pkls)

sstds = xr.open_dataset(SSTS)

xgr, ygr = np.mgrid[1:365:365j, -25:25:100j]
yinter = np.linspace(-25 - 0.5 / 2, 25 + 0.5 / 2, 101)
print(xgr[:, 0], ygr[0, ygr.shape[1]//2:], ygr[0, :ygr.shape[1]//2])

kdeobjs = [pk[('genday', 'genlat')][0] for pk in pkls]
[ko.set_bandwidth(0.18) for ko in kdeobjs]
nstms = [pk[('genday', 'genlat')][1] for pk in pkls]
print(kdeobjs[0])
print(nstms)

zvals = [ko.evaluate(np.vstack([xgr.ravel(), ygr.ravel()])).reshape(xgr.shape) * nstms[ii] for ii, ko in enumerate(kdeobjs)]
print(zvals[0])
znh = [zv[:, ygr.shape[1] // 2:] for zv in zvals] #must compute original kde with periodic boundary conditions
zsh = [np.roll(zv[:, :ygr.shape[1] // 2][:, ::-1], 182, axis=0) for zv in zvals]
zplt = [znh[ii] + zsh[ii] for ii in range(len(znh))]

#plt.contourf(xgr[:, ygr.shape[1] // 2:], ygr[:, ygr.shape[1] // 2:], zplt[1])
plt.rc('font', size=16)
plt.contourf(xgr, ygr, zvals[1])
plt.contour(xgr, ygr, zvals[1], levels=np.arange(.05, .4, .05))
plt.scatter(sstds.dayofyear, sstds['tos'].isel(case=1).idxmax('yh'), c='orange', marker='x')
plt.xticks(TICKDOY, [dt.strftime('%m-%d') for dt in TICKDATES], rotation=45)
#plt.show()
plt.close()

plt.rcParams['figure.figsize'] = (25, 5)
fig, axes = plt.subplots(1, 3, sharey=True, sharex=True)

for ii, ax in enumerate(axes):
   tocf = zvals[ii] if ii == 1 else zvals[ii] - zvals[1]
   cflvls = np.arange(.05, .46, .05) if ii == 1 else np.arange(-0.5, 0.51, .05)
   cmap = 'viridis' if ii == 1 else 'bwr'
   csf = ax.contourf(xgr, ygr, tocf, cmap=cmap, levels=cflvls)
   plt.colorbar(csf)
   ax.contour(xgr, ygr, zvals[1], colors='black', levels=np.arange(.05, .61, .05))
   ax.scatter(sstds.dayofyear, sstds['tos'].isel(case=ii).idxmax('yh'), c='aqua', marker='.')
   ax.set_xticks(TICKDOY, [dt.strftime('%m-%d') for dt in TICKDATES], rotation=45)
   ax.tick_params(axis='both', labelleft=True, right=True, top=True)
   ax.set_title(TTLS[ii][1])
   ax.set_title(TTLS[ii][0], loc='left')
   ax.set_xlabel('mm-dd')
   ax.set_ylabel('Latitude')

fig.tight_layout()
plt.savefig('TC_genlat_kde_sst.svg')
plt.show()

#dfs = [pd.read_parquet(pq) for pq in PQTS]
#print(dfs[0])
