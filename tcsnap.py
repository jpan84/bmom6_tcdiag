import os
import numpy as np
import xarray as xr
import cftime
import datetime
import matplotlib.pyplot as plt

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0606dlyout/'
HISTS = 'atm/hist/'
H1 = r'*h1.[0-9]*.nc'

pltdate = cftime.DatetimeNoLeap(50, 11, 9, hour=6)
LONBNDS = (55, 70)
LATBNDS = (15, 30)

#pltdate = cftime.DatetimeNoLeap(51, 5, 19, hour=18)
#LONBNDS = (233, 259)
#LATBNDS = (-35, -5)

def main():
   print('Opening history files...')
   ds = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H1)).sel(time=pltdate, lon=slice(*LONBNDS), lat=slice(*LATBNDS))

   contourfkwargs = {'cmap': 'YlGn', 'levels': np.arange(0, 121, 10), 'extend': 'max'}
   clabelkwargs = {'inline': 1, 'fontsize': 12, 'colors': 'black', 'fmt': '%d'}
   conlevs = np.arange(920, 1041, 4)

   cs = plt.contourf(ds.lon, ds.lat, ds.TMQ.values, **contourfkwargs)
   plt.colorbar(cs, ax=plt.gca())
   cs1 = plt.contour(ds.lon, ds.lat, ds.PSL.values / 100, levels=conlevs, colors='black')
   plt.gca().clabel(cs1, **clabelkwargs)
   plt.title('%s TMQ, PSL' % str(pltdate))
   plt.show()

   #print(dss[0].time.values)

   print('proj3.py done.')

#copied from https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
#weight monthly avgs to correctly compute annual avgs
def weighted_temporal_mean(ds, var):
   """
   weight by days in each month
   """
   # Determine the month length
   month_length = ds.time.dt.days_in_month

   # Calculate the weights
   wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()

   # Make sure the weights in each year add up to 1
   np.testing.assert_allclose(wgts.groupby("time.year").sum(xr.ALL_DIMS), 1.0)

   # Subset our dataset for our variable
   obs = ds[var]

   # Setup our masking for nan values
   cond = obs.isnull()
   ones = xr.where(cond, 0.0, 1.0)

   # Calculate the numerator
   obs_sum = (obs * wgts).resample(time="AS").sum(dim="time")

   # Calculate the denominator
   ones_out = (ones * wgts).resample(time="AS").sum(dim="time")

   # Return the weighted average
   return obs_sum / ones_out

if __name__ == '__main__':
   main()
