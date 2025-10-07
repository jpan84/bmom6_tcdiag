from sznl_funcs import stack_hemi_sznl
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


DIRI = './linevslat_h1i_0012-0014/'
TTL = ['UNSEED', 'CTRL', 'SEED']
FLDNM = 'TC_R4'

LNCLR = ['blue', 'orange']
YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.array([-90, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)

def main():
   if not os.path.exists(DIRI):
      os.makedirs(DIRI)

   ds = xr.open_dataset(os.path.join(DIRI, 'sznlzm.nc'))

   plt.rc('font', size=18)
   plt.rcParams['figure.figsize'] = (30, 6)
   subplot_kw = dict(xlim=(-1, 1))
   fig, axes = plt.subplot_mosaic([['(a)', '(b)', '(c)']], layout='constrained', sharey=True, subplot_kw=subplot_kw)

   for ii, (lbl, ax) in enumerate(axes.items()):
      [ax.axvline(yl, linewidth=0.5, color='gray') for yl in YLOC]
      ax.set_xlabel('Latitude [°]')
      ax.set_ylabel('%s fractional area' % FLDNM)
      toplt = ds[FLDNM].isel(case=ii)
      for jj, szn in enumerate(toplt.season):
         ax.plot(YSCL(toplt['latitudes']), toplt.sel(season=szn), color=LNCLR[jj], label='%s' % str(szn.data))
      ax.set_xticks(YLOC, YLAB)
      ax.set_title(TTL[ii])
      ax.set_title(lbl, loc='left')
      ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=12))

   plt.savefig(os.path.join(DIRI, 'replot_%s.png' % FLDNM))
   plt.close()

if __name__ == '__main__':
   main()
