import numpy as np
import os
import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/')
from paths import ARCHRT, ALIA, CTLIX, IXHORS
import xarray as xr
from sznl_funcs import stack_hemi_sznl, monthly2sznl
import matplotlib.pyplot as plt

TOTFIL = 'atm/uxzm_h0a_noncons_-90_90_0.25_TAUX.nc'
TCSFIL = 'atm/uxzm_h1i_noncons_-55_55.1_0.5_TC_R4_TAUX_PRECT.nc'
ALAT = 'latitudes'

TTLS = [ALIA[ii] + '$–$CTL' if ii != CTLIX else 'CTL' for ii in range(len(ALIA))]

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-60, 61, 10) #np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)

totdss = [xr.decode_cf(xr.open_dataset(os.path.join(ar, TOTFIL))) for ar in ARCHRT]
tcsdss = [xr.decode_cf(xr.open_dataset(os.path.join(ar, TCSFIL))) for ar in ARCHRT]

def main(pltvar='TAUX', pltsgn=-1):
   print([ds.time for ds in totdss])
   def agg_time(listdas, latnm='latitudes', diff=True):
      ymonmean = [da.groupby('time.month').mean() for da in listdas]
      twoszns = [stack_hemi_sznl(monthly2sznl(ym), antisym=False, latnm=latnm) for ym in ymonmean]
      halfyr = [ts.mean(dim='season') for ts in twoszns]
      if diff:
         return [hy - halfyr[CTLIX] if ii != CTLIX else hy for ii, hy in enumerate(halfyr)]
      return halfyr

   totplt = agg_time([pltsgn * ds[pltvar] for ds in totdss], latnm=ALAT)
   tcsplt = agg_time([pltsgn * ds[pltvar] for ds in tcsdss], latnm=ALAT)
   tcsplt[CTLIX] *= 5


   plt.rc('font', size=16)
   fig, axes = plt.subplots(2, 3, figsize=(22, 9), sharex=True)

   for ii, ax in enumerate(axes.ravel()):
       if ii == 4: continue
       ixh = IXHORS[ii]
       ax.plot(YSCL(totplt[ixh][ALAT]), totplt[ixh], label='mean')
       ax.plot(YSCL(tcsplt[ixh][ALAT]), tcsplt[ixh], label='TC R4')

       ax.set_title(TTLS[ixh])
       ax.set_title('(%s)' % chr(ord('a') + ii), loc='left')
       ax.axhline(0, c='gray')
       [ax.axvline(yl, c='gray', lw=0.5) for yl in YLOC]
       ax.set_xticks(YLOC, YLAB)
       ax.tick_params(top=True, right=True, labelbottom=True, labelleft=True)
       plt.legend()

       if not ixh == CTLIX:
          ax.set_xlim(YSCL(tcsplt[ixh][ALAT][0]), YSCL(tcsplt[ixh][ALAT][-1]))
          ax.set_ylim(-.015, .015)

   # Target the extra axis
   extra_ax = axes[1][1]
   extra_ax.set_axis_off() # Completely turn off ticks, labels, and borders safely

   # Create a clean solid fill box over the empty space
   extra_ax.add_patch(plt.Rectangle((0, 0.7), 1, 0.3, 
                                      facecolor='#000080', 
                                      transform=extra_ax.transAxes, 
                                      zorder=-1))

   # Add your centered text annotation
   extra_ax.text(
         0.5, 0.85, '↑↑ The orange line (TC R4) has been ↑↑\nmultiplied by 5 in panel (b) only.\nAll other lines show true values.',
         horizontalalignment='center',
         verticalalignment='center',
         transform=extra_ax.transAxes,
         fontsize=18,
         weight='bold', c='white'
   )

   fig.tight_layout()
   plt.savefig(f'TC_masked_{pltvar}.svg', bbox_inches='tight')
   plt.show()
   plt.close()


if __name__ == '__main__':
   main()
