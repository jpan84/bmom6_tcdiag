import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import os 
from paths import ALIA
from sznl_funcs import stack_hemi_sznl

FILI = '~/aquaptc/bmom6_tcdiag/paper1_post/TC_stats/260623_density_gmd/tcdens.nc'
SZN = 'JJASON'
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
fszn = lambda da, sz: da.sel(season=sz) if sz in da.season else da.sum('season')
shs = lambda da: stack_hemi_sznl(da, **stack_kw)
yvar = lambda da, case: da.sel(run=case) - int(case != 'ctrl') * da.sel(run='ctrl')
lysplt = fszn(shs(ds['lys']), SZN)
genplt = fszn(shs(ds['gen']), SZN)
aceplt = fszn(shs(ds['ace']), SZN)
fixplt = fszn(shs(ds['h6all']), SZN)
us_lys = fszn(shs(ds['us_lys']), SZN)
sd_gen = fszn(shs(ds['sd_gen']), SZN)
#events = [(-1 if TYPE[ii] == 'us' else 1) * fszn(shs(ds[str(cs.item()) + TYPE[ii]]), SZN) if TYPE[ii] is not None else 0 for ii, cs in enumerate(ds['run'])]
events = [fszn(shs(ds[str(cs.item()) + TYPE[ii]]), SZN) if TYPE[ii] is not None else 0 for ii, cs in enumerate(ds['run'])]
#print(events)

plt.rcParams['figure.figsize'] = (10, 22)
plt.rc('font', size=14)
fig, axes = plt.subplots(5, 2, sharex=True, sharey='col', subplot_kw=dict(xlim=(-1, 1)))
[ax.tick_params(labelbottom=True, labelleft=True, top=True, right=True) for ax in axes.ravel()]

for ii, ax in enumerate(axes[:, 0]):
   ax.plot(sinlat, yvar(lysplt, ds['run'][ii]), label='lysis' if TYPE[ii] is None else None)
   ax.plot(sinlat, yvar(genplt, ds['run'][ii]), label='genesis' if TYPE[ii] is None else None)
   if TYPE[ii] is None:
      dsin100mil = 1e8 / 2 / np.pi / (6.371e3)**2 #d(sinlat) to span 100 million sq km
      ax.annotate('', xy=(0, 7.), xytext=(dsin100mil, 7.), arrowprops=dict(arrowstyle=f'-[, widthB=0, lengthB=0',\
                      lw=1.5, color='purple', connectionstyle=f'bar,fraction={1/6}'), annotation_clip=False)
      ax.text(dsin100mil / 2, 9, '$100 \\times 10^6$ km$^2$', ha='center', va='bottom', fontsize=10, fontweight='bold', transform=ax.transData, c='purple')

   if TYPE[ii] == 'us':
      ax.plot(sinlat, us_lys.isel(run=ii), linestyle='dotted', color='C0', label='lysis at unseeding')
      #ax.plot(sinlat, yvar(lysplt, ds['run'][ii]) - us_lys.isel(run=ii), linestyle='dashed') #change in natural lysis points
      ax.plot(sinlat, events[ii], color='maroon', linestyle='dashdot', label='unseed events')

   if TYPE[ii] == 'sd':
      ax.plot(sinlat, sd_gen.isel(run=ii), linestyle='dotted', color='C1', label='seeded genesis')
      ax.plot(sinlat, yvar(genplt, ds['run'][ii]) - sd_gen.isel(run=ii), linestyle='dashed', color='C1', label='natural genesis') #change in natural gen points
      ax.plot(sinlat, events[ii], color='maroon', linestyle='dashdot', label='seed events')

   ax.set_yscale('symlog', linthresh=3.)
   ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
   ax.axhline(0, color='gray', linewidth=0.5)
   [ax.axvline(ll, color='gray', linewidth=0.5) for ll in LATLOC]
   [ax.axhline(hl, color='gray', linewidth=0.5) for hl in np.arange(-2, 7)]
   ax.set_yticks(np.concatenate((np.arange(-2.5, 2.5, .5), [3, 4, 6])))
   ax.set_ylim(-2.7, 6)
   ax.set_xticks(LATLOC2, LATLAB2)
   ax.set_title('(%s)' % chr(97 + 2 * ii), loc='left')

   ax.set_ylabel(ALIA[ii] + (' – CTL' if TYPE[ii] is not None else ''), fontsize=18)
   ax.legend(loc='best', fontsize=12)

for ii, ax in enumerate(axes[:, 1]):
   ax.plot(sinlat, yvar(aceplt, ds['run'][ii]), label='ACE [10$^4$ kt$^2$]', color='darkgreen')
   ax.plot(sinlat, yvar(fixplt, ds['run'][ii]), label='6-hourly TC fixes', color='black')

   ax.set_yscale('symlog', linthresh=30.)
   ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
   ax.axhline(0, color='gray', linewidth=0.5)
   [ax.axvline(ll, color='gray', linewidth=0.5) for ll in LATLOC]
   [ax.axhline(hl, color='gray', linewidth=0.5) for hl in np.arange(-20, 70, 10)]
   ax.set_yticks(np.concatenate((np.arange(-25, 30, 5), [30, 40, 60])))
   ax.set_ylim(-27, 60)
   ax.set_xticks(LATLOC2, LATLAB2)
   ax.set_title('(%s)' % chr(98 + 2 * ii), loc='left')

   ax.legend(loc='best', fontsize=12)


# --- Place this right after line 93 (after your loops finish) ---

from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox

# 1. Force a draw so Matplotlib calculates final subplot geometry
fig.canvas.draw() 

# 2. Extract the specific row axes (Row Index 2)
row_2_axes = [axes[2, 0], axes[2, 1]]

# 3. Combine their bounding boxes in figure coordinates
bboxes = [ax.get_position() for ax in row_2_axes]
row_bbox = Bbox.union(bboxes)

# 4. Add slight padding so the box doesn't tightly crowd the axis labels
padding = 0.02
padded_bbox = Bbox([[row_bbox.x0 - padding, row_bbox.y0 - padding],
                    [row_bbox.x1 + padding, row_bbox.y1 + padding]])

# 5. Create and style the unsaturated red bounding box
rect = Rectangle(
    (padded_bbox.x0, padded_bbox.y0), 
    padded_bbox.width, 
    padded_bbox.height,
    facecolor='#f8d7da',       # Unsaturated red / light pastel red
    edgecolor='#dc3545',       # Darker, muted red outline
    linewidth=1.5,
    linestyle='-',             # Solid border (or '--' for dashed)
    transform=fig.transFigure,  # Apply relative to the entire figure layout
    zorder=-1                   # Push to the absolute background layer
)

# 6. Make row 2 axes transparent so the background color peeks through
for ax in row_2_axes:
    ax.set_facecolor('none')

# 7. Add the bounding box to the figure canvas
fig.patches.append(rect)


fig.suptitle('TC density plots [%s$^{-1}$ (10$^6$ km$^2$)$^{-1}$]' % SZN)

fig.tight_layout()
plt.savefig('masterplot_%s.svg' % SZN)
plt.show()
plt.close()
