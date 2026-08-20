import numpy as np
import os
import sys
sys.path.insert(0, '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/')
sys.path.insert(0, '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post')
from paths import ARCHRT, ALIA, CTLIX, CASENAMES
import consts as c
import xarray as xr
from sznl_funcs import stack_hemi_sznl, monthly2sznl
import matplotlib.pyplot as plt

TOTP = '/glade/campaign/univ/upsu0032/jpan_aquaptc/%s/atm/uxzm_hist_h1i_noncons_-90.0_90.0_2.0_PRECT_TMQ.nc'
TCSP = '/glade/derecho/scratch/jpan/archive/nff_output/%s/uxzm_nff4mps_h1i_noncons_-90.0_90.0_2.0_TC_R4_PRECT.nc'
TOTE = '/glade/derecho/scratch/jpan/archive/nff_output/%s/hflso_subset_zonal_mean.nc'
TCSE = '/glade/derecho/scratch/jpan/archive/nff_output/%s/hflso_masked_4mps_zonal_mean.nc'

ALAT = 'latitudes'
OLAT = 'yh'

TTLS = [ALIA[ii] + '$–$CTL' if ii != CTLIX else 'CTL' for ii in range(len(ALIA))]

m2mm = 1e3
d2s = 8.64e4

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-60, 61, 10) #np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)

def main():
   def agg_time(listdas, latnm='latitudes', diff=True):
      ymonmean = [da.groupby('time.month').mean() for da in listdas]
      twoszns = [stack_hemi_sznl(monthly2sznl(ym), antisym=False, latnm=latnm) for ym in ymonmean]
      halfyr = [ts.mean(dim='season') for ts in twoszns]
      if diff:
         return [hy - halfyr[CTLIX] if ii != CTLIX else hy for ii, hy in enumerate(halfyr)]
      return halfyr

   totpds = [xr.open_dataset(TOTP % cs) for cs in CASENAMES]
   tcspds = [xr.open_dataset(TCSP % cs) for cs in CASENAMES]
   toteds = [xr.open_dataset(TOTE % cs).rename({OLAT: ALAT}) for cs in CASENAMES]
   tcseds = [xr.open_dataset(TCSE % cs).rename({OLAT: ALAT, 'hflso_masked': 'hflso'}) for cs in CASENAMES]

   totp_plt = agg_time([ds['PRECT'] * m2mm * d2s for ds in totpds], latnm=ALAT, diff=False)
   tcsp_plt = agg_time([ds['PRECT'] * m2mm * d2s for ds in tcspds], latnm=ALAT, diff=False)
   tcr4_plt = agg_time([ds['TC_R4'] for ds in tcspds], latnm=ALAT, diff=False)
   cwv_plt = agg_time([ds['TMQ'] for ds in totpds], latnm=ALAT, diff=False)

   tote_plt = agg_time([-ds['hflso'] / c.lv / c.rho_w * m2mm * d2s for ds in toteds], latnm=ALAT, diff=False)
   tcse_plt = agg_time([-ds['hflso'] / c.lv / c.rho_w * m2mm * d2s for ds in tcseds], latnm=ALAT, diff=False)
   tote_plt = [da.interp(latitudes=totp_plt[0][ALAT]) da in tote_plt]
   tcse_plt = [da.interp(latitudes=totp_plt[0][ALAT]) da in tcse_plt]

   plt.rcParams['figure.figsize'] = (8, 10)
   fig, axes = plt.subplots(3, 2, sharex=True)
   colors = plt.colormaps['inferno'](np.linspace(0.2, 0.8, 5))

   for ii in range(len(totp_plt)):
      axes[0][0].plot(totp_plt[ii][ALAT], totp_plt[ii], color=colors[ii])
      axes[0][0].plot(tcsp_plt[ii][ALAT], tcsp_plt[ii], color=colors[ii], ls='dashed')

      axes[0][1].plot(totp_plt[ii][ALAT], tcsp_plt[ii] / totp_plt[ii], color=colors[ii], ls='dashed')
      axes[0][1].plot(tcr4_plt[ii][ALAT], tcr4_plt[ii], color=colors[ii], ls='dotted')

      axes[1][0].plot(tote_plt[ii][ALAT], tote_plt[ii], color=colors[ii])
      axes[1][0].plot(tcse_plt[ii][ALAT], tcse_plt[ii], color=colors[ii], ls='dashed')

      axes[1][1].plot(tote_plt[ii][ALAT], tcse_plt[ii] / tote_plt[ii], color=colors[ii], ls='dashed')
      axes[1][1].plot(tcr4_plt[ii][ALAT], tcr4_plt[ii], color=colors[ii], ls='dotted')

   plt.show()

   exit()







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
       ax.set_xlabel('Latitude [°]')
       ax.set_ylabel('$\\tau_x$' + ('' if ii == 1 else ' anomaly') + ' [N m$^{-2}$]')
       ax.legend(loc='upper left')

       if not ixh == CTLIX:
          ax.set_xlim(YSCL(tcsplt[ixh][ALAT][0]), YSCL(tcsplt[ixh][ALAT][-1]))
          #ax.set_ylim(-.015, .015)

   # Target the extra axis
   extra_ax = axes[1][1]
   extra_ax.set_axis_off() # Completely turn off ticks, labels, and borders safely

   # Create a clean solid fill box over the empty space
   extra_ax.add_patch(plt.Rectangle((-0.1, 0.7), 1.2, 0.3, 
                                      facecolor='#000080', 
                                      transform=extra_ax.transAxes, 
                                      zorder=-1))

   # Add your centered text annotation
   extra_ax.text(
         0.5, 0.85, '↑ The orange line (TC R4) has been ↑\nmultiplied by 5 in panel (b) only.\nAll other lines show true values.',
         horizontalalignment='center',
         verticalalignment='center',
         transform=extra_ax.transAxes,
         fontsize=18,
         weight='bold', c='white'
   )

   fig.tight_layout()
   plt.savefig(f'TC_masked_{pltvar}.pdf', bbox_inches='tight')
   plt.show()
   plt.close()


if __name__ == '__main__':
   main()
