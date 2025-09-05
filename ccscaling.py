import os
import xarray as xr
import matplotlib.pyplot as plt

DIRI = 'linevslat_new_sznl_diff/'

lv = 2.501e6
Rv = 461.5

ds = xr.open_dataset(os.path.join(DIRI, 'sznlzm.nc'))
TS = ds['TS'].sel(season='SON')
QR = ds['TMQ'].sel(season='SON')
Tsq = TS.sel(case='CTRL') ** 2

#plt.plot(ds['latitudes'], lv / Rv / Tsq)
#plt.show()

ccfrac = lv / Rv / Tsq
dq_pred = TS.sel(case='UNSEED') * ccfrac * QR.sel(case='CTRL')
dq_true = QR.sel(case='UNSEED')

plt.plot(ds['latitudes'], dq_true, label='true')
plt.plot(ds['latitudes'], dq_pred, label='CC')
plt.title('$\Delta$TMQ UNSEED')
plt.axhline(0)
plt.legend()
plt.show()
