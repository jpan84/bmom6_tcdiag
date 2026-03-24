import pickle
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

DIRI = './0302_test_mf_vars/'
#FILI = 'UMF500_lat_UBOT_yy09dd01_05-35warm_allprec.nc'
#FILI = 'UMF500_SSTr_MSE850_yy09dd0x_05-35warm.nc'
MFIL = 'UMF500_SSTr_MSE850_yy09mm01-06_05-35warm.nc'
AFIL = 'areasr_SSTr_MSE850_yy09mm01-06_05-35warm.nc'
MFIL = 'UMF500_SSTr_MSE850r_yy06-08_05-35warm.nc'
AFIL = 'areasr_SSTr_MSE850r_yy06-08_05-35warm.nc'
#DIRO = 'mf_histo_btbc'
MVAR, XVAR, YVAR = MFIL.split('_')[:3]
AVAR, _, _ = AFIL.split('_')[:3]

TTLS = ['UNSEED$-$CTRL', 'CTRL', 'MSEED$-$CTRL']

mds = xr.open_dataset(os.path.join(DIRI, MFIL)).sum(dim='time')
ads = xr.open_dataset(os.path.join(DIRI, AFIL)).sum(dim='time')

###
#This code block is for ratio of fractions
###

mfrac = mds[MVAR] / mds[MVAR].sum()
afrac = ads[AVAR] / ads[AVAR].sum()
relcon = mfrac / afrac
relcon.loc[dict(case=['UNSEED', 'MSEED'])] -= relcon.sel(case='CTRL')

#SSTr, MSE850
subplot_kw = dict(ylim=(3.2, 3.8), xlim=(-5, 5))
subplot_kw = dict(ylim=(-4, 2.5), xlim=(-8, 5))
#pckw = [dict(cmap='seismic', vmin=-3e6, vmax=3e6), dict(cmap='nipy_spectral', vmin=0, vmax=3.5e7), dict(cmap='seismic', vmin=-3e6, vmax=3e6)]
pckw = [dict(cmap='nipy_spectral', norm=colors.LogNorm(vmin=1e-1, vmax=1e3)) for _ in range(3)] #raw
pckw = [dict(cmap='bwr', norm=colors.SymLogNorm(1e-2, vmin=-1e2, vmax=1e2)), dict(cmap='nipy_spectral', norm=colors.LogNorm(vmin=1e-1, vmax=1e2))] #diff
pckw.append(pckw[0])

plt.rc('font', size=14)
plt.rcParams['figure.figsize'] = [16, 5]

#deep tropical (Moist Tropics) case
fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, subplot_kw=subplot_kw, layout='constrained')
fig.suptitle('JJASON ratio of 500 hPa UMF to surface area [fraction / fraction]')
for ii, ax in enumerate(axes):
   pltda = relcon.isel(case=ii).T
   pc = ax.pcolormesh(mds[XVAR], mds[YVAR] / 1e4, pltda, shading='nearest', **pckw[ii])
   cb = plt.colorbar(pc, ax=ax, extend='both')

   ax.tick_params(right=True, top=True, labelleft=True)
   #ax.set_yticks(np.arange(3.2, 3.9, .1))
   ax.set_xlabel('SST\' [°C]')
   ax.set_ylabel('MSE850\' [10$^4$ J kg$^{-1}$]')
   ax.set_title(f'{TTLS[ii]}')
   ax.set_title(['(a)', '(b)', '(c)'][ii], loc='left')

#plt.show()
plt.close()

###
#This code block is for UMF normalized by area
###
relcon = mds[MVAR] / ads[AVAR]
relcon.loc[dict(case=['UNSEED', 'MSEED'])] -= relcon.sel(case='CTRL')

#SSTr, MSE850
subplot_kw = dict(ylim=(3.2, 3.8), xlim=(-5, 5))
subplot_kw = dict(ylim=(-4, 2.5), xlim=(-8, 5))
#pckw = [dict(cmap='seismic', vmin=-3e6, vmax=3e6), dict(cmap='nipy_spectral', vmin=0, vmax=3.5e7), dict(cmap='seismic', vmin=-3e6, vmax=3e6)]
pckw = [dict(cmap='nipy_spectral', norm=colors.LogNorm(vmin=1e10, vmax=1e13)) for _ in range(3)] #raw
pckw = [dict(cmap='bwr', norm=colors.SymLogNorm(1e9, vmin=-1e12, vmax=1e12)), dict(cmap='nipy_spectral', norm=colors.LogNorm(vmin=1e10, vmax=1e13))] #diff
pckw.append(pckw[0])

plt.rc('font', size=14)
plt.rcParams['figure.figsize'] = [16, 5]

#deep tropical (Moist Tropics) case
fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, subplot_kw=subplot_kw, layout='constrained')
fig.suptitle('JJASON ratio of 500 hPa UMF to surface area [fraction / fraction]')
for ii, ax in enumerate(axes):
   pltda = relcon.isel(case=ii).T
   pc = ax.pcolormesh(mds[XVAR], mds[YVAR] / 1e4, pltda, shading='nearest', **pckw[ii])
   cb = plt.colorbar(pc, ax=ax, extend='both')

   ax.tick_params(right=True, top=True, labelleft=True)
   #ax.set_yticks(np.arange(3.2, 3.9, .1))
   ax.set_xlabel('SST\' [°C]')
   ax.set_ylabel('MSE850\' [10$^4$ J kg$^{-1}$]')
   ax.set_title(f'{TTLS[ii]}')
   ax.set_title(['(a)', '(b)', '(c)'][ii], loc='left')

plt.show()
