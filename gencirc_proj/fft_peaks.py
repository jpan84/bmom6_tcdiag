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

peaks, _ = scipy.signal.find_peaks(ds[IXVAR], height=1.0, distance=3 * 4)
#print(peaks)
plt.plot(ds[IXVAR])
plt.scatter(peaks, ds[IXVAR].isel(time=peaks), c='orange')
plt.show()
