import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

DIRI = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr/atm/'
TFIL = '0012-0013_JJASON_onpres_1deg_zm1.nc'
FFIL = '0012-0013_JJASON_onpres_1deg_EPF_reorder.nc'
TVAR = 'T'
FVAR = 'EPz_EHF'
LATB = (10, 30)
PREB = (400, 700)

tds = xr.open_dataset(os.path.join(DIRI, TFIL)).squeeze()
fds = xr.open_dataset(os.path.join(DIRI, FFIL))

#print(tds.lat.values)
#print(tds.plev.values)
#print(fds.plev.values)

dTlat = tds[TVAR].interp(lat=LATB[1]) - tds[TVAR].interp(lat=LATB[0])
dTpav = (dTlat * dTlat['plev']).sum(dim='plev') / dTlat['plev'].sum()
#print(dTpav.mean())

fsel = fds[FVAR].sel(lat=slice(*LATB), plev=slice(*PREB))
latw = np.cos(np.deg2rad(fsel['lat']))
flav = (fsel * latw).sum(dim='lat') / latw.sum(dim='lat')
fpav = (flav * flav['plev']).sum(dim='plev') / flav['plev'].sum()

plt.rc('font', size=16)
plt.scatter(dTpav, fpav)
plt.title('Avg over 700–400 hPa')
plt.xlabel('$\Delta T_{30-10°N}$')
plt.ylabel('$F^{(p)}$ [kg s-2]')
plt.axhline(0, c='gray')
plt.savefig('EPz_vs_dTy.png', bbox_inches='tight')
plt.show()
