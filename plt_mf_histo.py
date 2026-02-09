import pickle
import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

FILI = 'UMF500_0012dd01_07-20warm.nc'
#FILI = 'umf500_0012dd01_20-33warm.nc'
DIRO = 'mf_histo_btbc'
TTLS = ['UNSEED$-$CTRL', 'CTRL', 'SEED$-$CTRL']

ds = xr.open_dataset(FILI).sum(dim='time') #forgot to change time normalization in 1st script when retaining time dim so need to sum here instead of avg
umfdif = ds['UMF500'].copy()
umfdif.loc[dict(case=['UNSEED', 'SEED'])] -= umfdif.sel(case='CTRL')
print(umfdif.max())
umfdif /= ds['xwidth'].data[None, :, None]
umfdif /= ds['ywidth'].data[None, None, :]


#RAWVLIM = dict(vmin=0, vmax=1.5e11)
#DIFSZNVLIM = dict(vmin=-1.5e11, vmax=1.5e11)
#DIFCASEVLIM = dict(vmin=-5e10, vmax=5e10)

plt.rc('font', size=14)
plt.rcParams['figure.figsize'] = [16, 5]

#deep tropical (Moist Tropics) case
fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, subplot_kw=dict(xlim=(3.2, 3.8), ylim=(29, 37)), layout='constrained')
fig.suptitle('JJASON 500 hPa UMF [kg s$^{{-1}}$ (°C J kg$^{{-1}}$)$^{-1}$]')
for ii, ax in enumerate(axes):
   rawda = ds['UMF500'].isel(case=ii)
   pckw = dict(cmap='nipy_spectral', vmin=0, vmax=3.5e7) if ii == 1 else dict(cmap='seismic', vmin=-2e7, vmax=2e7) if ii == 2\
             else dict(cmap='seismic', vmin=-4e6, vmax=4e6)
   pltda = umfdif.isel(case=ii).T
   pc = ax.pcolormesh(ds['MSE850'] / 1e5, ds['SST'] - 273.15, pltda, shading='nearest', **pckw)
   cb = plt.colorbar(pc, ax=ax, extend='both')
   #if ii != 1:
   #   cb.ax.set_yscale('log')

   ax.tick_params(right=True, top=True, labelleft=True)
   ax.set_xlabel('MSE 850 [$10^5$ J kg$^{-1}$]')
   ax.set_ylabel('SST [°C]')
   ax.set_title(f'{TTLS[ii]}\nsum: {rawda.sum():.2e} kg s$^{{-1}}$')
   ax.set_title(['(a)', '(b)', '(c)'][ii], loc='left')

plt.show()

#subtropical (Subtropical Transition Zone) case
fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, subplot_kw=dict(xlim=(3.2, 3.8), ylim=(26, 34)), layout='constrained')
fig.suptitle('JJASON 500 hPa UMF [kg s$^{{-1}}$ (°C J kg$^{{-1}}$)$^{-1}$]')
for ii, ax in enumerate(axes):
   rawda = ds['UMF500'].isel(case=ii)
   pckw = dict(cmap='nipy_spectral', vmin=0, vmax=3.5e7) if ii == 1 else dict(cmap='seismic', vmin=-2e7, vmax=2e7) if ii == 2\
             else dict(cmap='seismic', vmin=-4e6, vmax=4e6)
   pltda = umfdif.isel(case=ii).T
   pc = ax.pcolormesh(ds['MSE850'] / 1e5, ds['SST'] - 273.15, pltda, shading='nearest', **pckw)
   cb = plt.colorbar(pc, ax=ax, extend='both')
   #if ii != 1:
   #   cb.ax.set_yscale('log')

   ax.tick_params(right=True, top=True, labelleft=True)
   ax.set_xlabel('MSE 850 [$10^5$ J kg$^{-1}$]')
   ax.set_ylabel('SST [°C]')
   ax.set_title(f'{TTLS[ii]}\nsum: {rawda.sum():.2e} kg s$^{{-1}}$')
   ax.set_title(['(d)', '(e)', '(f)'][ii], loc='left')

plt.show()

exit()
#plot every raw histogram
for ii, al in enumerate(ds['case']):
   #for jj, hm in enumerate(HEM):
   rawda = ds['UMF500'].sel(case=al)
   pltda = umfdif.sel(case=al)
   plt.pcolormesh(ds['MSE850'] / 1e5, ds['SST'] - 273.15, pltda, shading='nearest', cmap='viridis', **RAWVLIM)
   plt.xlim(*XLIM)
   plt.ylim(*YLIM)
   plt.xlabel('MSE 850 [J kg$^{-1}$]')
   plt.ylabel('SST [°C]')
   plt.title(f'{al} 500 hPa UMF [kg s$^{{-1}}$]\n{rawda.sum():.3e}')
   plt.colorbar(extend='max')

   plt.show()
      #plt.savefig(os.path.join(DIRO, 'umf_binned_2d_%s_%s.png' % (al, hm)), bbox_inches='tight')
      #plt.close()

#plot every seasonal difference
for ii, al in enumerate(ALIS):
   dif = histos[ii * 2][0].T - histos[ii * 2 + 1][0].T
   plt.pcolormesh(MSE850, SST, dif, shading='flat', cmap='seismic', **DIFSZNVLIM)
   plt.xlim(*XLIM)
   plt.ylim(*YLIM)
   plt.xlabel('MSE 850 [J kg$^{-1}$]')
   plt.ylabel('SST [°C]')
   plt.title(f'Warm minus cool 500 hPa UMF [kg s$^{{-1}}$]\n{dif.sum():.3e}')
   plt.colorbar(extend='both')
   plt.savefig(os.path.join(DIRO, 'umf_binned_2d_%s_warm-cool.png' % al), bbox_inches='tight')
   plt.close()

#plot every run minus ctrl
for ii, al in enumerate(ALIS):
   for jj, hm in enumerate(HEM):
      dif = histos[ii * 2 + jj][0].T - histos[2 + jj][0].T
      plt.pcolormesh(MSE850, SST, dif, shading='flat', cmap='seismic', **DIFCASEVLIM)
      plt.xlim(*XLIM)
      plt.ylim(*YLIM)
      plt.xlabel('MSE 850 [J kg$^{-1}$]')
      plt.ylabel('SST [°C]')
      plt.title(f'{al} minus CTRL 500 hPa UMF [kg s$^{{-1}}$]\n{dif.sum():.3e}')
      plt.colorbar(extend='both')
      plt.savefig(os.path.join(DIRO, 'umf_binned_2d_%s-CTRL_%s.png' % (al, hm)), bbox_inches='tight')
      plt.close()
