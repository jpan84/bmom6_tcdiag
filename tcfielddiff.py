import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

CTLDIR = './tcfields4mps_250417_ctrl/'
EXPDIR = './tcfields4mps_250415_unseed/'
FLD = ['hflso_neg', 'hflso_pos']
SGN = -1

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-50, 51, 10)

ctlall = xr.open_dataset(os.path.join(CTLDIR, 'means_all.nc'))
ctltcs = xr.open_dataset(os.path.join(CTLDIR, 'means_tcs.nc'))
expall = xr.open_dataset(os.path.join(EXPDIR, 'means_all.nc'))
exptcs = xr.open_dataset(os.path.join(EXPDIR, 'means_tcs.nc'))
dss = [ctlall, expall, ctltcs, exptcs]

flds = [SGN * sum([ds[dv] for dv in FLD]) for ds in dss]
diff = [flds[1] - flds[0], flds[-1] - flds[-2]]

plt.rc('font', size=16)
sinlat = YSCL(ctlall['yh'])
for ii, szn in enumerate(diff[0]['season']):
   plt.plot(sinlat, diff[0].sel(season=szn), color=['blue', 'orange'][ii], label=str(szn) + '_all')
   plt.plot(sinlat, diff[1].sel(season=szn), color=['blue', 'orange'][ii], linestyle='dashed', label=str(szn) + '_TCs')
   plt.xlim(YSCL(YLAB[0]), YSCL(YLAB[-1]))
   plt.xticks(YSCL(YLAB), YLAB)
   plt.axhline(0, linestyle='dashed', color='black')

   plt.ylim(-8, 3.5)
   plt.yticks(np.arange(-8, 4))
   plt.title('LHFLX')

plt.show()
