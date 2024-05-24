import os
import numpy as np
import xarray as xr
import cftime
import datetime
import matplotlib.pyplot as plt
import scipy.stats as stats
tdel = datetime.timedelta

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0515ctrl/'
#ARCHV = '/glade/derecho/scratch/youweima/archive/'
#CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.yr500_branch/'
#setsid nohup qcmd -q casper -A UPSU0064 -l walltime=00:30:00  -l select=1:ncpus=10:mem=200GB  python3 proj3.py
HISTS = 'atm/hist/'
H0 = r'*h0.[0-9]*.nc'

def main():
   print('Opening history files...')
   ds = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H0))
   ds = ds.assign_coords(coords=dict(time=ds.time - tdel(days=1))) #timestamp of e.g., 0001-01.nc is Feb 1, so subtract a month from time coord
   #print(dss[0].time.values)

   assert(isinstance(ds.time.values[0], cftime._cftime.DatetimeNoLeap))

   print('Computing annual means...')
   ts_ann = weighted_temporal_mean(ds, 'TS')
   restom = weighted_temporal_mean(ds, 'FSNT') - weighted_temporal_mean(ds, 'FLNT')

   #print(ts_ann[0].shape)

   print('Computing global means...')
   latwgt = np.repeat(ds.gw.isel(time=0).values[:, np.newaxis], ds.lon.shape[0], axis=-1)
   print(latwgt.shape)
   globav = lambda da: (latwgt * da).sum(dim=['lat', 'lon']) / latwgt.sum()
   ts_av = globav(ts_ann)
   rad_av = globav(restom)

   '''
   print('Computing sfc temp anoms...')
   #tbase = ts_av[0].sel(time=slice(ts_ann[0].time[0].values + tdel(days = 365 * 20),)).mean(dim='time')
   tbase = ts_av[0][dict(time=slice(20,50))].mean(dim='time')
   print(ts_av[0][dict(time=slice(20,50))].time.values)
   #print(ts_ann[0].time[0] + tdel(days = 365 * 20))
   tanoms = [ts - tbase for ts in ts_av]
   '''

   print('Plotting time series...')
   plt.rc('font', size=16)
   plt.plot(ts_av)
   #plt.hlines(0, 0, 50, linestyle='--', color='black')
   plt.xlabel('Year')
   plt.ylabel('$T_s$ [K]')
   #plt.legend()
   plt.savefig('ts.png', bbox_inches='tight')
   plt.close()

   plt.plot(rad_av)
   plt.hlines(0, 0, 20, linestyle='--', color='black')
   plt.xlabel('Year')
   plt.ylabel('Top of model\nnet radiative forcing [W m$^{-2}$]')
   plt.legend()
   plt.savefig('restom.png', bbox_inches='tight')
   plt.close()

   '''
   print('Making Gregory plots...')
   plt.scatter(tanoms[1], rad_av[1], c='orange')
   plt.xlabel('$\Delta T_s$ [K]')
   plt.ylabel('Top of model\nnet radiative forcing [W m$^{-2}$]')
   plt.title('%s capelmt = 150' % LBLS[1])
   plt.savefig('gregory.png', bbox_inches='tight')

   print('Computing ECS...')
   linreg = stats.linregress(rad_av[1], tanoms[1])
   print('Slope: %f' % linreg.slope)
   print('Intercept: %f' % linreg.intercept)
   rr = np.arange(0, 10, 1)
   tt = linreg.intercept + linreg.slope * rr
   plt.plot(tt, rr, linestyle='--', color='black')
   plt.text(10, 8, '2ECS (intercept): %.2f\nSlope: %.2f' % (linreg.intercept, linreg.slope))
   plt.savefig('gregory.png', bbox_inches='tight')
   '''

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
