import pickle
import os
import matplotlib.pyplot as plt

FILI = 'umf500_0012dd0x_histo_mse850_SST.pkl'
DIRO = 'mf_histo_SST'

RAWVLIM = dict(vmin=0, vmax=1.5e11)
DIFSZNVLIM = dict(vmin=-1.5e11, vmax=1.5e11)
DIFCASEVLIM = dict(vmin=-5e10, vmax=5e10)

XLIM = (3.0e5, 3.8e5)
YLIM = (27, 37.1)

if not os.path.exists(DIRO):
   os.makedirs(DIRO)

histos = None
with open(FILI, 'rb') as fl:
   histos = pickle.load(fl)

#print(histos)
ALIS = histos.pop(0)
HEM = histos.pop(0)

MSE850 = histos[0][1]
SST = histos[0][2] - 273.15

plt.rc('font', size=14)

#plot every raw histogram
for ii, al in enumerate(ALIS):
   for jj, hm in enumerate(HEM):
      plt.pcolormesh(MSE850, SST, histos[ii * 2 + jj][0].T, shading='flat', cmap='viridis', **RAWVLIM)
      plt.xlim(*XLIM)
      plt.ylim(*YLIM)
      plt.xlabel('MSE 850 [J kg$^{-1}$]')
      plt.ylabel('SST [°C]')
      plt.title(f'{al} 500 hPa UMF [kg s$^{{-1}}$]\n{histos[ii * 2 + jj][0].sum():.3e}')
      plt.colorbar(extend='max')
      plt.savefig(os.path.join(DIRO, 'umf_binned_2d_%s_%s.png' % (al, hm)), bbox_inches='tight')
      plt.close()

#plot every seasonal difference
for ii, al in enumerate(ALIS):
   dif = histos[ii * 2][0].T - histos[ii * 2 + 1][0].T
   plt.pcolormesh(MSE850, SST, dif, shading='flat', cmap='bwr', **DIFSZNVLIM)
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
      plt.pcolormesh(MSE850, SST, dif, shading='flat', cmap='bwr', **DIFCASEVLIM)
      plt.xlim(*XLIM)
      plt.ylim(*YLIM)
      plt.xlabel('MSE 850 [J kg$^{-1}$]')
      plt.ylabel('SST [°C]')
      plt.title(f'{al} minus CTRL 500 hPa UMF [kg s$^{{-1}}$]\n{dif.sum():.3e}')
      plt.colorbar(extend='both')
      plt.savefig(os.path.join(DIRO, 'umf_binned_2d_%s-CTRL_%s.png' % (al, hm)), bbox_inches='tight')
      plt.close()
