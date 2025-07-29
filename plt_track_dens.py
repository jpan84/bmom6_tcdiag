import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from sznl_funcs import stack_hemi_sznl

FILI = '/glade/u/home/jpan/aquaptc/tempest/250725_density_sznl/tcdens.nc'

LATLAB = np.arange(-90, 90.1, 15).astype(np.int_)
LATLAB = np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90])
LATLOC = np.sin(np.deg2rad(LATLAB))

ds = xr.open_dataset(FILI)
sinlat = np.sin(np.deg2rad(ds['lat']))
print(ds)

plt.rcParams['figure.figsize'] = (20, 6)
plt.rc('font', size=16)
fig, axes = plt.subplots(2, 3, sharex=True, sharey='row', subplot_kw=dict(xlim=(-1, 1)))

stack_kw = dict(antisym=False, sznnm='season', latnm='lat')
lysplt = stack_hemi_sznl(ds['lys'], **stack_kw).sum(dim='season')
genplt = stack_hemi_sznl(ds['gen'], **stack_kw).sum(dim='season')
aceplt = stack_hemi_sznl(ds['ace'], **stack_kw).sum(dim='season')
allplt = stack_hemi_sznl(ds['h6all'], **stack_kw).sum(dim='season')

axes[0][0].set_title('UNSEED')
axes[0][0].plot(sinlat, lysplt.isel(run=0) - lysplt.isel(run=1))
axes[0][0].plot(sinlat, genplt.isel(run=0) - genplt.isel(run=1))
axes[0][0].set_ylim(-3.5, 3.5)
axes[0][0].set_xticks(LATLOC, LATLAB)

axes[0][1].set_title('CTRL')
axes[0][1].plot(sinlat, lysplt.isel(run=1), label='lysis')
axes[0][1].plot(sinlat, genplt.isel(run=1), label='genesis')
axes[0][1].legend(loc='best')

axes[0][2].set_title('SEED')
axes[0][2].plot(sinlat, lysplt.isel(run=2) - lysplt.isel(run=1))
axes[0][2].plot(sinlat, genplt.isel(run=2) - genplt.isel(run=1))


axes[1][0].plot(sinlat, aceplt.isel(run=0) - aceplt.isel(run=1))
axes[1][0].plot(sinlat, allplt.isel(run=0) - allplt.isel(run=1), linestyle='dotted', color='black')
axes[1][0].set_ylim(-25, 65)
axes[1][0].set_xticks(LATLOC, LATLAB)

axes[1][1].plot(sinlat, aceplt.isel(run=1), label='ACE [10$^4$ kt$^2$]')
axes[1][1].plot(sinlat, allplt.isel(run=1), label='6-hourly TC fixes', linestyle='dotted', color='black')
axes[1][1].legend(loc='best')
axes[1][1].set_xlabel('Latitude')

axes[1][2].plot(sinlat, aceplt.isel(run=2) - aceplt.isel(run=1))
axes[1][2].plot(sinlat, allplt.isel(run=2) - allplt.isel(run=1), linestyle='dotted', color='black')

fig.suptitle('TC density plots [yr$^{-1}$ (10$^6$ km$^2$)$^{-1}$]\nJJASON')

fig.tight_layout()
plt.show()
