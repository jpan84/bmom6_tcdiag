import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

VAR = 'EHF500' #main var of interest
LATNM = 'latitudes'

DIRI = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/EHF_EMF_500_warmNH_0012-0014_eofs/'
EIGS = '251117_NH_eigvals.nc'
EOFS = '251117_NH_eofs.nc'
MEAN = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/EHF_EMF_500_warmNH_0012-0014_timmean.nc'

eids = xr.open_dataset(os.path.join(DIRI, EIGS)).squeeze()
eods = xr.open_dataset(os.path.join(DIRI, EOFS)).squeeze()
meds = xr.open_dataset(MEAN).squeeze()

totvar = eids[VAR].sum(dim='time')
expvar = eids[VAR] / totvar

plt.rc('font', size=16)
plt.scatter(range(1, expvar.size + 1), expvar)
plt.xlim(.8, 6.2)
[plt.axhline(fr, lw=0.5, c='gray') for fr in np.arange(0, 0.7, .1)]
plt.xlabel('EOF mode')
plt.ylabel('Variance explained')
plt.savefig('expvar.png', bbox_inches='tight')
#plt.show()
plt.close()

plt.plot(meds[LATNM], meds[VAR], c='black', lw=2.5)
[plt.axvline(ll, lw=0.5, c='gray') for ll in np.arange(0, 30, 5)]
plt.ylabel('Mean EHF [K m/s]')
plt.xlabel('Lat [°]')
ax1 = plt.gca().twinx()
[plt.plot(eods[LATNM], eods[VAR].isel(time=ii), label=ii + 1) for ii in range(4)]
plt.ylabel('EOF patterns')
#plt.show()
plt.legend()
plt.savefig('EOF500_patterns.png', bbox_inches='tight')
plt.close()
