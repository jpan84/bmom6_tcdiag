import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import os 
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sznl_funcs import stack_hemi_sznl

FILI = '/glade/u/home/jpan/aquaptc/tempest/260415_density_5exp/tcdens.nc'
SZN = 'JJA'
TYPE = ['us', 'us', None, 'sd', 'sd']

YSCL = lambda deglat: np.sin(np.deg2rad(deglat))
LATLAB = np.arange(-90, 90.1, 15).astype(np.int_)
LATLAB = np.array([-90, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 90])
LATLOC = np.sin(np.deg2rad(LATLAB))
LATLAB2 = np.arange(-90, 91, 30)
LATLOC2 = YSCL(LATLAB2)

MT = YSCL(np.array([10, 20]))
STZ = YSCL(np.array([20, 30]))

ds = xr.open_dataset(FILI)
sinlat = np.sin(np.deg2rad(ds['lat']))
print(ds)

stack_kw = dict(antisym=False, sznnm='season', latnm='lat')
shs = lambda da: stack_hemi_sznl(da, **stack_kw).sel(season=SZN)
yvar = lambda da, case: da.sel(run=case) - int(case != 'ctrl') * da.sel(run='ctrl')
lysplt = shs(ds['lys'])
genplt = shs(ds['gen'])
aceplt = shs(ds['ace'])
fixplt = shs(ds['h6all'])
us_lys = shs(ds['us_lys'])
sd_gen = shs(ds['sd_gen'])
events = [(-1 if TYPE[ii] == 'us' else 1) * shs(ds[str(cs.item()) + TYPE[ii]]) if TYPE[ii] is not None else 0 for ii, cs in enumerate(ds['run'])]
events = [shs(ds[str(cs.item()) + TYPE[ii]]) if TYPE[ii] is not None else 0 for ii, cs in enumerate(ds['run'])]
#print(events)

plt.rcParams['figure.figsize'] = (10, 20)
plt.rc('font', size=14)
fig, axes = plt.subplots(5, 2, sharex=True, sharey='col', subplot_kw=dict(xlim=(-1, 1)))
[ax.tick_params(labelbottom=True, labelleft=True, right=True) for ax in axes.ravel()]

for ii, ax in enumerate(axes[:, 0]):
   ax.plot(sinlat, yvar(lysplt, ds['run'][ii]), label='lysis' if TYPE[ii] is None else None)
   ax.plot(sinlat, yvar(genplt, ds['run'][ii]), label='genesis' if TYPE[ii] is None else None)

   if TYPE[ii] == 'us':
      ax.plot(sinlat, us_lys.isel(run=ii), linestyle='dotted', color='C0')
      #ax.plot(sinlat, yvar(lysplt, ds['run'][ii]) - us_lys.isel(run=ii), linestyle='dashed') #change in natural lysis points
      ax.plot(sinlat, events[ii], color='maroon', linestyle='dashdot', label='unseed events')

   if TYPE[ii] == 'sd':
      ax.plot(sinlat, sd_gen.isel(run=ii), linestyle='dotted', color='C1')
      ax.plot(sinlat, yvar(genplt, ds['run'][ii]) - sd_gen.isel(run=ii), linestyle='dashed', color='C1') #change in natural gen points
      ax.plot(sinlat, events[ii], color='maroon', linestyle='dashdot', label='seed events')

   ax.set_yscale('symlog', linthresh=1.2)
   ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
   ax.axhline(0, color='gray', linewidth=0.5)
   [ax.axvline(ll, color='gray', linewidth=0.5) for ll in LATLOC]
   ax.set_yticks(np.concatenate((np.arange(-1, 1.5, .25), np.arange(2, 4))))
   ax.set_ylim(-1.1, 3.5)
   ax.set_xticks(LATLOC2, LATLAB2)
   ax.set_title('(%s)' % chr(97 + 2 * ii), loc='left')

   ax.set_ylabel(str(ds['run'][ii].item()), fontsize=18)
   ax.legend(loc='best', fontsize=12)

for ii, ax in enumerate(axes[:, 1]):
   ax.plot(sinlat, yvar(aceplt, ds['run'][ii]), label='ACE [10$^4$ kt$^2$]' if TYPE[ii] is None else None, color='darkgreen')
   ax.plot(sinlat, yvar(fixplt, ds['run'][ii]), label='6-hourly TC fixes' if TYPE[ii] is None else None, color='black')

   ax.set_yscale('symlog', linthresh=12)
   ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
   ax.axhline(0, color='gray', linewidth=0.5)
   [ax.axvline(ll, color='gray', linewidth=0.5) for ll in LATLOC]
   ax.set_yticks(np.concatenate((np.arange(-10, 13, 2), np.arange(15, 35, 5))))
   ax.set_ylim(-12, 31)
   ax.set_xticks(LATLOC2, LATLAB2)
   ax.set_title('(%s)' % chr(98 + 2 * ii), loc='left')

   ax.legend(loc='best', fontsize=12)

fig.suptitle('%s TC density plots [yr$^{-1}$ (10$^6$ km$^2$)$^{-1}$]' % SZN)

fig.tight_layout()
plt.savefig('/glade/u/home/jpan/aquaptc/tempest/260415_density_5exp/masterplot_%s.png' % SZN)

exit()

plt.rcParams['figure.figsize'] = (16, 8)
plt.rc('font', size=14)
fig, axes = plt.subplots(2, 3, sharex=True, sharey='row', subplot_kw=dict(xlim=(-1, 1)))
[ax.tick_params(labelbottom=True, labelleft=True, right=True) for ax in axes.ravel()]

stack_kw = dict(antisym=False, sznnm='season', latnm='lat')
lysplt = stack_hemi_sznl(ds['lys'], **stack_kw).sel(season=SZN)#.sum(dim='season')
genplt = stack_hemi_sznl(ds['gen'], **stack_kw).sel(season=SZN)#.sum(dim='season')
aceplt = stack_hemi_sznl(ds['ace'], **stack_kw).sel(season=SZN)#.sum(dim='season')
unsplt = -stack_hemi_sznl(ds['unseeds'], **stack_kw).sel(season=SZN)
sdsplt = stack_hemi_sznl(ds['seeds'], **stack_kw).sel(season=SZN)
allplt=stack_hemi_sznl(ds['h6all'], **stack_kw).sel(season=SZN)#.sum(dim='season')

axes[0][0].set_title('UNSEED–CTRL')
axes[0][0].plot(sinlat, lysplt.isel(run=0) - lysplt.isel(run=1))
axes[0][0].plot(sinlat, genplt.isel(run=0) - genplt.isel(run=1))
axes[0][0].plot(sinlat, unsplt, color='maroon', linestyle='dashdot', label='unseed events')
axes[0][0].axhline(0, color='gray', linewidth=0.5)
[axes[0][0].axvline(ll, color='gray', linewidth=0.5) for ll in LATLOC]
axes[0][0].set_yticks(np.arange(-1, 1.1, .25))
axes[0][0].set_ylim(-1.3, 1.3)
axes[0][0].set_xticks(LATLOC2, LATLAB2)
axes[0][0].set_title('(a)', loc='left')
axes[0][0].legend(loc='best', fontsize=12)

axes[0][1].set_title('CTRL')
axes[0][1].plot(sinlat, lysplt.isel(run=1), label='lysis')
axes[0][1].plot(sinlat, genplt.isel(run=1), label='genesis')
axes[0][1].axhline(0, color='gray', linewidth=0.5)
[axes[0][1].axvline(ll, color='gray', linewidth=0.5) for ll in LATLOC]
axes[0][1].legend(loc='best', fontsize=12)
axes[0][1].set_title('(b)', loc='left')
#axes[0][1].annotate('', xy=(MT[0], 1.35), xytext=(MT[1], 1.35), arrowprops=dict(arrowstyle=f'-[, widthB=0, lengthB=0',\
#                        lw=1.5, color='purple', connectionstyle=f'bar,fraction={1/6}'), annotation_clip=False)
#axes[0][1].text(sum(MT) / 2, 1.6, 'MT', ha='center', va='bottom', fontsize=10, fontweight='bold', transform=axes[0][1].transData, c='purple')
#axes[0][1].annotate('', xy=(STZ[0], 1.35), xytext=(STZ[1], 1.35), arrowprops=dict(arrowstyle=f'-[, widthB=0, lengthB=0',\
#                        lw=1.5, color='#E4D00A', connectionstyle=f'bar,fraction={1/6}'), annotation_clip=False)
#axes[0][1].text(sum(STZ) / 2, 1.6, 'STZ', ha='center', va='bottom', fontsize=10, fontweight='bold', transform=axes[0][1].transData, c='#E4D00A')

axes[0][2].set_title('MSEED–CTRL')
axes[0][2].plot(sinlat, lysplt.isel(run=2) - lysplt.isel(run=1))
axes[0][2].plot(sinlat, genplt.isel(run=2) - genplt.isel(run=1))
axes[0][2].plot(sinlat, sdsplt, color='maroon', linestyle='dashdot', label='seed events')
axes[0][2].axhline(0, color='gray', linewidth=0.5)
[axes[0][2].axvline(ll, color='gray', linewidth=0.5) for ll in LATLOC]
axes[0][2].set_title('(c)', loc='left')
axes[0][2].legend(loc='best', fontsize=12)

axes[1][0].plot(sinlat, aceplt.isel(run=0) - aceplt.isel(run=1), color='darkgreen')
axes[1][0].plot(sinlat, allplt.isel(run=0) - allplt.isel(run=1), linestyle='dashed', color='black')

axes[1][0].set_xticks(LATLOC2, LATLAB2)
axes[1][0].axhline(0, color='gray', linewidth=0.5)
[axes[1][0].axvline(ll, color='gray', linewidth=0.5) for ll in LATLOC]
axes[1][0].set_title('(d)', loc='left')

axes[1][1].plot(sinlat, aceplt.isel(run=1), label='ACE [10$^4$ kt$^2$]', color='darkgreen')
axes[1][1].plot(sinlat, allplt.isel(run=1), label='6-hourly TC fixes', linestyle='dashed', color='black')
axes[1][1].legend(loc='best', fontsize=12)
axes[1][1].set_xlabel('Latitude')
axes[1][1].axhline(0, color='gray', linewidth=0.5)
[axes[1][1].axvline(ll, color='gray', linewidth=0.5) for ll in LATLOC]
axes[1][1].set_title('(e)', loc='left')

axes[1][2].plot(sinlat, aceplt.isel(run=2) - aceplt.isel(run=1), color='darkgreen')
axes[1][2].plot(sinlat, allplt.isel(run=2) - allplt.isel(run=1), linestyle='dashed', color='black')
axes[1][2].axhline(0, color='gray', linewidth=0.5)
[axes[1][2].axvline(ll, color='gray', linewidth=0.5) for ll in LATLOC]
axes[1][2].set_title('(f)', loc='left')

[ax.axvspan(*MT, fill=True, fc='purple', zorder=2, alpha=0.08) for ax in axes.ravel()]
[ax.axvspan(*STZ, fill=True, fc='yellow', zorder=2, alpha=0.12) for ax in axes.ravel()]

fig.suptitle('%s TC density plots [yr$^{-1}$ (10$^6$ km$^2$)$^{-1}$]' % SZN)

fig.tight_layout()
plt.savefig('/glade/u/home/jpan/aquaptc/tempest/260324_density/masterplot_%s.png' % SZN)
plt.show()
plt.close()
