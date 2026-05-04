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
SSTS = 'sstmaxlat_ydaymean.nc'
TTLS = [('(a)', 'UNSEED$-$CTRL'), ('(b)', 'CTRL'), ('(c)', 'SEED$-$CTRL')]
NYR = 10. #maybe update to be more precise (TODO don't hardcode)

TICKDATES = [cftime.datetime(1, ii, 1, calendar='noleap') for ii in range(1, 13)]
TICKDOY = [int(dt.strftime('%-j')) for dt in TICKDATES]

YSCL = lambda y: np.sin(np.deg2rad(y))

pkls = []
for fl in FILI:
   with open(os.path.join(DIRI, fl), 'rb') as ofl:
      pkls.append(pickle.load(ofl))
#print(pkls)

sstds = xr.open_dataset(SSTS)

latlim = 25.
#xgr, ygr = np.mgrid[1:365:365j, -25:25:100j]
#yinter = np.linspace(-25 - 0.5 / 2, 25 + 0.5 / 2, 101)
yinter = np.linspace(-YSCL(latlim + 0.5), YSCL(latlim + 0.5), 101)
xlin = np.linspace(1., 365, 365)
dday = xlin[1] - xlin[0]
yarea = 2 * np.pi * 6.371e6**2 / 1e6 / 1e3**2 * np.diff(yinter) #million sq km
ymid = (yinter[:-1] + yinter[1:]) / 2
xgr, ygr = np.meshgrid(xlin, ymid, indexing='ij')
yplt = np.rad2deg(np.arcsin(ygr))
print(xgr[:, 0], ygr[0, ygr.shape[1]//2:], ygr[0, :ygr.shape[1]//2])

kdeobjs = [pk[('genday', 'genmu')][0] for pk in pkls]
[ko.set_bandwidth(0.08) for ko in kdeobjs]
nstms = [pk[('genday', 'genmu')][1] for pk in pkls]
print(kdeobjs[0])
print(nstms)

zvals = [ko.evaluate(np.vstack([xgr.ravel(), ygr.ravel()])).reshape(xgr.shape) * nstms[ii] / NYR / yarea[None, :] for ii, ko in enumerate(kdeobjs)] #stms/yr/doy/million sq km
raw_integs = [(zv * dday * yarea[None, :]).sum() for zv in zvals]
renorm = [nstms[ii] / NYR / ri for ii, ri in enumerate(raw_integs)]
zvals = [zv * renorm[ii] for ii, zv in enumerate(zvals)]
print(zvals[0])
znh = [zv[:, ygr.shape[1] // 2:] for zv in zvals] #must compute original kde with periodic boundary conditions
zsh = [np.roll(zv[:, :ygr.shape[1] // 2][:, ::-1], 182, axis=0) for zv in zvals]
zplt = [znh[ii] + zsh[ii] for ii in range(len(znh))]
print(yarea.sum())
print((zvals[1] * dday * yarea[None, :]).sum())

#plt.contourf(xgr[:, ygr.shape[1] // 2:], ygr[:, ygr.shape[1] // 2:], zplt[1])
plt.rc('font', size=16)
plt.contourf(xgr, yplt, zvals[1])
plt.contour(xgr, yplt, zvals[1], levels=np.arange(.05, .4, .05))
#plt.scatter(sstds.dayofyear, sstds['tos'].isel(case=1).idxmax('yh'), c='orange', marker='x')
plt.xticks(TICKDOY, [dt.strftime('%m-%d') for dt in TICKDATES], rotation=45)
#plt.show()
plt.close()

plt.rcParams['figure.figsize'] = (25, 5)
fig, axes = plt.subplots(1, 3, sharey=True, sharex=True)

for ii, ax in enumerate(axes):
   tocf = zvals[ii] if ii == 1 else zvals[ii] - zvals[1]
   cflvls = np.arange(0, 2.1e-2, 2e-3) if ii == 1 else np.arange(-2e-2, 2.1e-2, 4e-3)
   cmap = 'viridis' if ii == 1 else 'bwr'
   csf = ax.contourf(xgr, yplt, tocf, cmap=cmap, levels=cflvls)
   plt.colorbar(csf)
   ax.contour(xgr, yplt, zvals[1], colors='black', levels=np.arange(4e-3, 2e-2, 4e-3))

   #tmaxlat = sstds['tos'].isel(case=ii).idxmax('yh')
   tmaxlat = sstds['tos'].isel(case=ii)
   print(tmaxlat.values)
   #exit()
   flipdoys = tmaxlat.dayofyear.isel(dayofyear=np.where(np.sign(tmaxlat.data[1:]) - np.sign(tmaxlat.data[:-1]))[0])
   ax.scatter(sstds.dayofyear, tmaxlat, c='aqua', marker='.')
   [ax.axvline(fd, c='gray', linestyle='dashed') for fd in flipdoys]
   fddt = [cftime.num2date(fd, units='days since 0000-12-31', calendar='noleap') for fd in flipdoys]
   fdmmdd = [ft.strftime('%m-%d') for ft in fddt]
   [ax.text(fd - 10, 0, fdmmdd[jj]) for jj, fd in enumerate(flipdoys)]

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
print(zvals[ii] * yarea.sum())
