import pickle
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

FILI = 'UMF500_0012dd01_07-20warm.nc'
#FILI = 'umf500_0012dd01_20-33warm.nc'
FILI = 'UMF850_OM850_OM500_0007dd1x_00-30warm.nc'
#FILI = 'UMF850_FLUT_OM500_0012dd0x_00-30warm.nc'
FILI = 'UMF850_OM850_OM500_yy0607_00-30warm.nc'
#FILI = 'UMF850_OM850_OM500_yy0607_00-30warm_allprec.nc'
DIRI = './0302_test_mf_vars/'
FILI = 'UMF500_lat_UBOT_yy09dd01_05-35warm_allprec.nc'
#DIRO = 'mf_histo_btbc'
THEVAR, XVAR, YVAR = FILI.split('_')[:3]
TTLS = ['UNSEED$-$CTRL', 'CTRL', 'MSEED$-$CTRL']

ds = xr.open_dataset(os.path.join(DIRI, FILI)).sum(dim='time') #forgot to change time normalization in 1st script when retaining time dim so need to sum here instead of avg
umfdif = ds[THEVAR].copy()
umfdif.loc[dict(case=['UNSEED', 'MSEED'])] -= umfdif.sel(case='CTRL')
print(umfdif.max())
if 'lat' in ds.dims:
   print('lat normalization')
   sinlat2 = np.sin(np.deg2rad(ds['lat'].data + ds['xwidth'].data / 2))
   sinlat1 = np.sin(np.deg2rad(ds['lat'].data - ds['xwidth'].data / 2))
   umfdif = umfdif / (2 * np.pi * 6.371e6**2 * (sinlat2 - sinlat1))[None, :, None]
else:
   umfdif /= ds['xwidth'].data[None, :, None]
umfdif /= ds['ywidth'].data[None, None, :]

#umfdif = -umfdif

#RAWVLIM = dict(vmin=0, vmax=1.5e11)
#DIFSZNVLIM = dict(vmin=-1.5e11, vmax=1.5e11)
#DIFCASEVLIM = dict(vmin=-5e10, vmax=5e10)
subplot_kw = dict(xlim=(3.2, 3.8), ylim=(29, 37))
pckw = [dict(cmap='seismic', vmin=-4e6, vmax=4e6), dict(cmap='nipy_spectral', vmin=0, vmax=3.5e7), dict(cmap='seismic', vmin=-2e7, vmax=2e7)]
subplot_kw = dict(xlim=(-1, 0), ylim=(-.5, 0.5))
pckw = [dict(cmap='seismic', vmin=-5e11, vmax=5e11), dict(cmap='nipy_spectral', vmin=0, vmax=2e13), dict(cmap='seismic', vmin=-5e11, vmax=5e11)]
#subplot_kw = dict(xlim=(-2, .5), ylim=(-2, .5))
#pckw = [dict(cmap='seismic', vmin=-2e11, vmax=2e11), dict(cmap='nipy_spectral', vmin=0, vmax=1.5e12), dict(cmap='seismic', vmin=-2e11, vmax=2e11)]

#subplot_kw = dict(xlim=(80, 340), ylim=(-.5, .5))
#pckw = [dict(cmap='seismic', vmin=-1.5e10, vmax=1.5e10), dict(cmap='nipy_spectral', vmin=0, vmax=7.5e10), dict(cmap='seismic', vmin=-1.5e10, vmax=1.5e10)]

#lat, SST
subplot_kw = dict(xlim=(5, 35), ylim=(26, 37))
pckw = [dict(cmap='seismic', vmin=-4e-3, vmax=4e-3), dict(cmap='nipy_spectral', vmin=0, vmax=1e-2), dict(cmap='seismic', vmin=-4e-3, vmax=4e-3)]

#lat, UBOT
subplot_kw = dict(xlim=(5, 35), ylim=(0, 40))
pckw = [dict(cmap='seismic', vmin=-5e-4, vmax=5e-4), dict(cmap='nipy_spectral', vmin=0, vmax=1e-3), dict(cmap='seismic', vmin=-5e-4, vmax=5e-4)]

plt.rc('font', size=14)
plt.rcParams['figure.figsize'] = [16, 5]

#deep tropical (Moist Tropics) case
fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, subplot_kw=subplot_kw, layout='constrained')
fig.suptitle('JJASON 500 hPa UMF [kg s$^{{-1}}$ (°C J kg$^{{-1}}$)$^{-1}$]')
fig.suptitle('JJASON 850 hPa UMF [kg s$^{{-1}}$ (Pa s$^{{-1}}$)$^{-2}$]')
#fig.suptitle('JJASON 850 hPa UMF [kg s$^{{-1}}$ (Pa s$^{{-1}}$ W m$^{{-2}}$)$^{-1}$]')
fig.suptitle('JJASON 500 hPa UMF [kg s$^{{-1}}$ (°C m s$^{{-1}}$)$^{-1}$]')
for ii, ax in enumerate(axes):
   rawda = ds[THEVAR].isel(case=ii)
   pltda = umfdif.isel(case=ii).T
   #pc = ax.pcolormesh(ds['MSE850'] / 1e5, ds['SST'] - 273.15, pltda, shading='nearest', **pckw[ii])
   #pc = ax.pcolormesh(ds[XVAR], ds[YVAR], pltda, shading='nearest', **pckw[ii])
   pc = ax.pcolormesh(ds[XVAR], ds[YVAR] - 273.15, pltda, shading='nearest', **pckw[ii])
   pc = ax.pcolormesh(ds[XVAR], ds[YVAR], pltda, shading='nearest', **pckw[ii])
   cb = plt.colorbar(pc, ax=ax, extend='both')
   #ax.plot([max(subplot_kw['xlim'][0], subplot_kw['ylim'][0]), min(subplot_kw['xlim'][1], subplot_kw['ylim'][1])], [max(subplot_kw['xlim'][0], subplot_kw['ylim'][0]), min(subplot_kw['xlim'][1], subplot_kw['ylim'][1])], c='green')
   #if ii != 1:
   #   cb.ax.set_yscale('log')

   ax.tick_params(right=True, top=True, labelleft=True)
   ax.set_xlabel('$\omega_{850}$ [Pa s$^{-1}$]')
   #ax.set_xlabel('OLR [W m$^{-2}$]')
   ax.set_xlabel('lat [deg]')
   ax.set_ylabel('$\omega_{500}$ [Pa s$^{-1}$]')
   ax.set_ylabel('SST [°C]')
   ax.set_ylabel('UBOT speed [m/s]')
   ax.set_title(f'{TTLS[ii]}\nsum: {rawda.sum():.2e} kg s$^{{-1}}$')
   ax.set_title(['(a)', '(b)', '(c)'][ii], loc='left')

plt.show()
