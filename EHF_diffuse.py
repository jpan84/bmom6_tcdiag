import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

DIRI = 'streamf_sznl/'
F1 = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl__TEM.nc'
F2 = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl_b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1_TEM.nc'

LATLAB = np.arange(-90, 91, 30)
LATLOC = np.sin(np.deg2rad(LATLAB))

ds1 = xr.open_dataset(os.path.join(DIRI, F1))
ds2 = xr.open_dataset(os.path.join(DIRI, F2))
print(ds1)

K_dif = -ds1['vpTHp'] / ds1['TH_y']
print(np.quantile(K_dif, np.arange(0, 1.01, .1)))
plt.contourf(ds1['lat'], ds1['plev'], K_dif.sel(season='JJA'), extend='both', levels=np.arange(-4e4, 2.1e5, 2e4), norm=colors.TwoSlopeNorm(0), cmap='PuOr_r')
plt.gca().invert_yaxis()
plt.colorbar()
plt.savefig(os.path.join(DIRI, 'Kdif_ctrl.png'))
plt.close()

plt.contourf(ds1['lat'], ds1['plev'], ds2['vpTHp'].sel(season='JJA'), extend='both', norm=colors.TwoSlopeNorm(0), cmap='bwr', levels=np.arange(-2, 2.1, 0.4))
plt.gca().invert_yaxis()
plt.colorbar()
#plt.show()
plt.close()

contourkwargs = {'colors': 'black', 'levels': 2.**np.arange(-2, 7, 1)}
contourkwargs['levels'] = np.concatenate((-contourkwargs['levels'][::-1], contourkwargs['levels']))
EHF_diff = -ds2['TH_y'] * K_dif
print(np.quantile(EHF_diff, np.arange(0, 1.01, .1)))
csf = plt.contourf(np.sin(np.deg2rad(ds1['lat'])), ds1['plev'], EHF_diff.sel(season='JJA'), extend='both', norm=colors.TwoSlopeNorm(0), cmap='bwr', levels=np.arange(-.2, .21, .04))
plt.contour(np.sin(np.deg2rad(ds1['lat'])), ds1['plev'], ds2['PSI_vT'].sel(season='JJA') / 1e9, **contourkwargs)
plt.xticks(LATLOC, LATLAB)
plt.yticks(np.arange(100, 1001, 100))
plt.ylim(1000, 100)
plt.yscale('log')
#plt.gca().invert_yaxis()
plt.colorbar(csf)
plt.savefig(os.path.join(DIRI, 'EHFdiff_Kdif.png'))
#plt.show()
plt.close()

K2 = -(ds1['vpTHp'] + ds2['vpTHp']) / (ds1['TH_y'] + ds2['TH_y'])
delK = K2 - K_dif
csf = plt.contourf(np.sin(np.deg2rad(ds1['lat'])), ds1['plev'], delK.sel(season='JJA'), extend='both', norm=colors.TwoSlopeNorm(0), cmap='PuOr_r', levels=np.arange(-5e4, 5.1e4, 1e4))
plt.contour(np.sin(np.deg2rad(ds1['lat'])), ds1['plev'], ds2['PSI_vT'].sel(season='JJA') / 1e9, **contourkwargs)
plt.xticks(LATLOC, LATLAB)
plt.yticks(np.arange(100, 1001, 100))
plt.ylim(1000, 100)
plt.yscale('log')
plt.colorbar(csf)
plt.savefig(os.path.join(DIRI, 'deltaK.png'))
#plt.show()
plt.close()

EHF_diff = -ds1['TH_y'] * delK
print(np.quantile(EHF_diff, np.arange(0, 1.01, .1)))
csf = plt.contourf(np.sin(np.deg2rad(ds1['lat'])), ds1['plev'], EHF_diff.sel(season='JJA'), extend='both', norm=colors.TwoSlopeNorm(0), cmap='bwr', levels=np.arange(-4, 4.1, .5))
plt.contour(np.sin(np.deg2rad(ds1['lat'])), ds1['plev'], ds2['PSI_vT'].sel(season='JJA') / 1e9, **contourkwargs)
plt.xticks(LATLOC, LATLAB)
plt.yticks(np.arange(100, 1001, 100))
plt.ylim(1000, 100)
plt.yscale('log')
#plt.gca().invert_yaxis()
plt.colorbar(csf)
plt.savefig(os.path.join(DIRI, 'EHFdiff_delK.png'))
#plt.show()
plt.close()
