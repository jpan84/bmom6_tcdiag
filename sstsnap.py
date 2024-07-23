import os
import numpy as np
import xarray as xr
import cftime
import datetime
import matplotlib.pyplot as plt

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0711dlyout/'
HISTS = 'ocn/hist/'
H1 = r'*hmd.[0-9]*.nc'
H1 = '*hmd*1030_11_*.nc'
H1 = r'*hmd*1029_0[7-8]_*.nc'
H1 = '*hmd*1031_02_*.nc'

var = 'oml'
#var = 'tos'

#strong TC
pltdate = cftime.DatetimeNoLeap(1030, 11, 20, hour=18)
LONBNDS = (50, 60)
LATBNDS = (15, 25)

#weak TC
pltdate = cftime.DatetimeNoLeap(1029, 8, 4, hour=6)
LONBNDS = (35, 45)
LATBNDS = (16, 26)

#strong ETC
pltdate = cftime.DatetimeNoLeap(1031, 2, 13, hour=0)
LONBNDS = (176.25, 186.25)
LATBNDS = (-57.3, -47.3)

#intensifying TC

TBNDS = (pltdate - datetime.timedelta(days=7), pltdate + datetime.timedelta(days=7))

#pltdate = cftime.DatetimeNoLeap(51, 5, 19, hour=18)
#LONBNDS = (233, 259)
#LATBNDS = (-35, -5)

def main():
   print('Opening history files...')
   ds = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H1)).sel(time=slice(*TBNDS), xh=slice(*LONBNDS), yh=slice(*LATBNDS))

   trc = ds[var].mean(dim=['xh', 'yh'])
   taxis = [(tt - pltdate).days + (tt - pltdate).seconds/86400 for tt in trc.time.values]

   print(taxis)

   plt.plot(taxis, trc.values)
   plt.title('%s %s' % (str(pltdate), var))
   plt.show()

if __name__ == '__main__':
   main()
