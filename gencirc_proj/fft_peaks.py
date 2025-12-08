import xarray as xr
import scipy
import numpy as np
import matplotlib.pyplot as plt

IXFIL = '0012-0013_JJASON_EEHF_ix_std.nc'
IXVAR = 'eehf_ix_std'

ds = xr.open_dataset(IXFIL)

print(ds['time'][ds.time.size // 2 - 3: ds.time.size // 2 + 3])
gb = ds.groupby(ds['time'].dt.year)

freqs = None
pwrs = []
fs_per_yr = 365. / 0.25
for yr, ds_yr in gb:
   print(yr)
   #print(scipy.fft.fft(ds_yr[IXVAR]))

   fq, psd = scipy.signal.periodogram(ds_yr[IXVAR], fs=fs_per_yr, scaling='density')
   freqs = fq
   pwrs.append(psd)
   #print(freqs[:5], psd[:5])

plt.rc('font', size=16)
pd_yr = 1. / freqs
pspec = np.mean(pwrs, axis=0)
plt.plot(fq, pspec)

pd_ticks = np.array([365, 90, 30, 14, 7, 1, 0.5]) #days
fr_ticks = 1 / (pd_ticks / 365) #cycle per year
print(pd_ticks)
print(fr_ticks)
#plt.gca().set_xticks(fr_ticks)
ax1 = plt.gca().twiny()
ax1.set_xticks(fr_ticks, labels=pd_ticks)
plt.show()
plt.close()

sepdays = 3
stdthresh = 1.0
peaks, _ = scipy.signal.find_peaks(ds[IXVAR], height=stdthresh, distance=sepdays * 4)
trofs, _ = scipy.signal.find_peaks(-ds[IXVAR], height=stdthresh, distance=sepdays * 4)
print(peaks.size, trofs.size)
pkda = xr.DataArray(ds['time'].isel(time=peaks).data, dims=['time1'], attrs=dict(sepdays=sepdays, stdthresh=stdthresh))
#pkda = pkda.assign_coords(time=pkda.data)
print(pkda)
tfda = xr.DataArray(ds['time'].isel(time=trofs).data, dims=['time2'], attrs=dict(sepdays=sepdays, stdthresh=stdthresh))
#tfda = tfda.assign_coords(time=tfda.data)
xr.Dataset(data_vars=dict(peaks=pkda, troughs=tfda)).to_netcdf('ehf_events_sep%d_sig%.1f.nc' % (sepdays, stdthresh))

plt.plot(ds[IXVAR])
plt.scatter(peaks, ds[IXVAR].isel(time=peaks), c='orange')
plt.scatter(trofs, ds[IXVAR].isel(time=trofs), c='green')
plt.show()
