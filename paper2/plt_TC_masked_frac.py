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
from matplotlib.lines import Line2D

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
   cwv_plt = agg_time([ds['TMQ'] / c.rho_w * m2mm for ds in totpds], latnm=ALAT, diff=False)

   tote_plt = agg_time([-ds['hflso'] / c.lv / c.rho_w * m2mm * d2s for ds in toteds], latnm=ALAT, diff=False)
   tcse_plt = agg_time([-ds['hflso'] / c.lv / c.rho_w * m2mm * d2s for ds in tcseds], latnm=ALAT, diff=False)
   tote_plt = [da.interp(latitudes=totp_plt[0][ALAT]) for da in tote_plt]
   tcse_plt = [da.interp(latitudes=totp_plt[0][ALAT]) for da in tcse_plt]

   plt.rcParams['figure.figsize'] = (8, 8)
   fig, axes = plt.subplots(3, 2, sharex=True)
   colors = plt.colormaps['cividis'](np.linspace(0.2, 0.8, 5))

   for ii in range(len(totp_plt)):
      axes[0][0].plot(YSCL(totp_plt[ii][ALAT]), totp_plt[ii], color=colors[ii], label=ALIA[ii])
      axes[0][0].plot(YSCL(tcsp_plt[ii][ALAT]), tcsp_plt[ii], color=colors[ii], ls='dashed')
      axes[0][0].set_ylabel('$P$ [mm d$^{-1}$]')
      axes[0][0].set_xticks(YLOC, YLAB)
      axes[0][0].set_xlim(YLOC[0], YLOC[-1])
      axes[0][0].legend()

      axes[0][1].plot(YSCL(totp_plt[ii][ALAT]), tcsp_plt[ii] / totp_plt[ii], color=colors[ii], ls='dashed')
      axes[0][1].plot(YSCL(tcr4_plt[ii][ALAT]), tcr4_plt[ii], color=colors[ii], ls='dotted')
      custom_handles = [
        Line2D([], [], color='gray', linestyle='--', label='$P$'),
        Line2D([], [], color='gray', linestyle=':', label='$r_4$')
      ]
      axes[0][1].set_ylabel('Fraction of zonal mean')
      axes[0][1].legend(handles=custom_handles)

      axes[1][0].plot(YSCL(tote_plt[ii][ALAT]), tote_plt[ii], color=colors[ii])
      axes[1][0].plot(YSCL(tcse_plt[ii][ALAT]), tcse_plt[ii], color=colors[ii], ls='dashed')
      axes[1][0].set_ylabel('$E$ [mm d$^{-1}$]')

      axes[1][1].plot(YSCL(tote_plt[ii][ALAT]), tcse_plt[ii] / tote_plt[ii], color=colors[ii], ls='dashed')
      axes[1][1].plot(YSCL(tcr4_plt[ii][ALAT]), tcr4_plt[ii], color=colors[ii], ls='dotted')

      axes[2][0].plot(YSCL(tote_plt[ii][ALAT]), totp_plt[ii] - tote_plt[ii], color=colors[ii])
      axes[2][0].plot(YSCL(tcse_plt[ii][ALAT]), tcsp_plt[ii] - tcse_plt[ii], color=colors[ii], ls='dashed')
      axes[2][0].set_ylabel('$P-E$ [mm d$^{-1}$]')

      axes[2][1].plot(YSCL(tote_plt[ii][ALAT]), (tcsp_plt[ii] - tcse_plt[ii]) / cwv_plt[ii], color=colors[ii], ls='dashed')
      axes[2][1].plot(YSCL(tote_plt[ii][ALAT]), (totp_plt[ii] - tote_plt[ii]) / cwv_plt[ii], color=colors[ii], ls='dotted')
      axes[2][1].set_ylabel('Normalized moisture sink $(P-E)/CWV$ [d$^{-1}$]')

   [ax.tick_params(top=True, labelbottom=True, right=True) for ax in axes.ravel()]
   [ax.set_title('(%s)' % chr(ord('a') + ii), loc='left') for ii, ax in enumerate(axes.ravel())]
   fig.tight_layout()
   plt.savefig('abs_P_E_frac.png', bbox_inches='tight')
   plt.close()

   totp_dif = [totp_plt[ii] - totp_plt[CTLIX] for ii in range(len(totp_plt))]
   tcsp_dif = [tcsp_plt[ii] - tcsp_plt[CTLIX] for ii in range(len(tcsp_plt))]
   tote_dif = [tote_plt[ii] - tote_plt[CTLIX] for ii in range(len(tote_plt))]
   tcse_dif = [tcse_plt[ii] - tcse_plt[CTLIX] for ii in range(len(tcse_plt))]
   totp_pme = [totp_dif[ii] - tote_dif[ii] for ii in range(len(totp_dif))]
   tcsp_pme = [tcsp_dif[ii] - tcse_dif[ii] for ii in range(len(tcsp_dif))]

   plt.rcParams['figure.figsize'] = (12, 6)
   fig, axes = plt.subplots(2, 3, sharex=True)

   for ii in range(len(totp_dif)):
      if ii == CTLIX: continue
      row = 0 if ii <= CTLIX else 1

      #axes[0][0].plot(YSCL(totp_plt[ii][ALAT]), totp_dif[ii], color=colors[ii], label=ALIA[ii])
      #axes[0][0].set_ylabel('$P$ [mm d$^{-1}$]')
      #axes[0][0].set_xticks(YLOC, YLAB)
      #axes[0][0].set_xlim(YLOC[0], YLOC[-1])
      #axes[0][0].legend()

      #axes[1][0].plot(YSCL(tcsp_plt[ii][ALAT]), tcsp_dif[ii], color=colors[ii], ls='dashed')

      axes[row][0].plot(YSCL(totp_plt[ii][ALAT]), totp_dif[ii], color=colors[ii], label=ALIA[ii])
      axes[row][0].plot(YSCL(tcsp_plt[ii][ALAT]), tcsp_dif[ii], color=colors[ii], ls='dashed')
      axes[row][0].set_ylabel('$P$ [mm d$^{-1}$]')
      axes[row][0].set_xticks(YLOC, YLAB)
      axes[row][0].set_xlim(YLOC[0], YLOC[-1])
      axes[row][0].legend()

      axes[row][1].plot(YSCL(tote_plt[ii][ALAT]), tote_dif[ii], color=colors[ii], label=ALIA[ii])
      axes[row][1].plot(YSCL(tcse_plt[ii][ALAT]), tcse_dif[ii], color=colors[ii], ls='dashed')
      axes[row][1].set_ylabel('$E$ [mm d$^{-1}$]')
      axes[row][1].set_xticks(YLOC, YLAB)
      axes[row][1].set_xlim(YLOC[0], YLOC[-1])
      axes[row][1].legend()

      axes[row][2].plot(YSCL(totp_pme[ii][ALAT]), totp_pme[ii], color=colors[ii], label=ALIA[ii])
      axes[row][2].plot(YSCL(tcsp_pme[ii][ALAT]), tcsp_pme[ii], color=colors[ii], ls='dashed')
      axes[row][2].set_ylabel('$P-E$ [mm d$^{-1}$]')
      axes[row][2].set_xticks(YLOC, YLAB)
      axes[row][2].set_xlim(YLOC[0], YLOC[-1])
      axes[row][2].legend()

   [ax.tick_params(top=True, labelbottom=True, right=True) for ax in axes.ravel()]
   [ax.set_title('(%s)' % chr(ord('a') + ii), loc='left') for ii, ax in enumerate(axes.ravel())]
   [ax.axhline(0, lw=0.5, c='gray') for ax in axes.ravel()]
   fig.tight_layout()
   plt.savefig('dif_P_E.png', bbox_inches='tight')
   plt.show()
   plt.close()



if __name__ == '__main__':
   main()
