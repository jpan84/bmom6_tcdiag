import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

MO_WTs = np.array([30, 31, 31, 30, 31, 30])
ticklevs = np.array([50, 70, 100, 200, 300, 500, 700, 850])
FILI = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr/atm/hist_onpres_1deg/ymonmean_VT_zm.nc'
FILI_h1 = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr/atm/0012-0013_JJASON_onpres_1deg_zm.nc'

ds = xr.open_mfdataset(FILI).squeeze()
ds = ds.sel(time=(ds['time'].dt.month.isin(np.arange(6, 12))))
ehf = ds['VT'] - ds['V'] * ds['T']
ehf_tm = (ehf * MO_WTs[:, None, None]).sum(dim='time') / MO_WTs.sum()

bbox = (ds['plev'] <= 8.5e4) & (ds['plev'] >= 3e4) & (ds['lat'] >= 0) & (ds['lat'] <= 35)
isneg = ehf_tm < 0

levels = 2. ** np.arange(-.5, 5.1, 0.5)
levels = np.concatenate((-levels[::-1], [0], levels))
plt.contourf(ds['lat'], ds['plev'], ehf_tm, cmap='RdBu_r', levels=levels, norm=colors.SymLogNorm(0.5))
plt.xlim(-10, 60)
plt.ylim(8.5e4, 5e3)
plt.yscale('log')
#plt.yticks(np.arange(1e4, 1.01e5, 1e4), np.arange(100, 1001, 100))
plt.yticks(ticklevs * 100, ticklevs)
plt.title('JJASON 10-yr mean $v\'T\'$ [K m/s]')
plt.xlabel('lat')
plt.ylabel('p [hPa]')
plt.gca().tick_params(top=True, right=True)
cb = plt.colorbar()
cb.ax.set_yticks(levels[::2])
plt.contour(ds['lat'], ds['plev'], bbox, color='black')
#plt.show()
plt.savefig('EHF_JJASON_10yr.png')
plt.close()

h1ds = xr.open_mfdataset(FILI_h1)
ehf_h1 = h1ds['VT'] - ds['V'] * ds['T']
masswt = ds['plev'] * np.cos(np.deg2rad(ds['lat']))
eehf_ix = (ehf_tm * ehf_h1 * masswt * bbox * isneg).sum(['plev', 'lat'])
stdix = (eehf_ix - eehf_ix.mean(dim='time')) / eehf_ix.std(dim='time')
with ProgressBar():
   xr.Dataset(data_vars=dict(eehf_ix_std=stdix)).to_netcdf('0012-0013_JJASON_EEHF_ix_std.nc')
