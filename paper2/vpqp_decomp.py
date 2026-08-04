import os
import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/')
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post')
from sznl_funcs import stack_hemi_sznl, monthly2sznl
from paths import ARCHRT, ALIA, CTLIX, IXHORS
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

#DIRI = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm'
totfil = 'atm/uxzm_hist_h1i_noncons_-90.0_90.0_2.0_TAUX_PRECT_V850_Q850_V850.Q850.nc'
tcsfil = 'atm/uxzm_nff4mps_h1i_noncons_-90.0_90.0_2.0_TAUX_PRECT_V850_Q850_V850.Q850.nc'
bkgfil = 'atm/uxzm_nff4mpsinvert_h1i_noncons_-90.0_90.0_2.0_TAUX_PRECT_V850_Q850_V850.Q850.nc'
ALAT = 'latitudes'

TTLS = [ALIA[ii] + '$–$CTL' if ii != CTLIX else 'CTL' for ii in range(len(ALIA))]

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-60, 61, 10) #np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)

def rey_dec_vq(vds, qds, vvar='V850', qvar='Q850'):
   vq = None
   if vds is qds:
      vq = vds[vvar + '_dot_' + qvar]
   else:
      vq = xr.zeros_like(vds[vvar]) #the pointwise product is zero for complementary masks (e.g., v_TC, q_BG)
   vq_mean = vds[vvar] * qds[qvar]
   vq_eddy = vq - vq_mean
   return [vq, vq_mean, vq_eddy]  

def agg_time(listdas, latnm='latitudes', antisym=False, diff=True):
   ymonmean = [da.groupby('time.month').mean() for da in listdas]
   twoszns = [stack_hemi_sznl(monthly2sznl(ym), antisym=antisym, latnm=latnm) for ym in ymonmean]
   halfyr = [ts.mean(dim='season') for ts in twoszns]
   if diff:
      return [hy - halfyr[CTLIX] if ii != CTLIX else hy for ii, hy in enumerate(halfyr)]
   return halfyr

def main():
   totdss = [xr.open_dataset(os.path.join(ar, totfil)) for ar in ARCHRT]
   tcsdss = [xr.open_dataset(os.path.join(ar, tcsfil)) for ar in ARCHRT]
   bkgdss = [xr.open_dataset(os.path.join(ar, bkgfil)) for ar in ARCHRT]

   covs = []
   for tods, tcds, bgds in zip(totdss, tcsdss, bkgdss):
      vq_tot = rey_dec_vq(tods, tods)[-1]
      vq_tcs = rey_dec_vq(tcds, tcds)[-1]
      vq_bkg = rey_dec_vq(bgds, bgds)[-1]
      vq_t2b = rey_dec_vq(tcds, bgds)[-1]
      vq_b2t = rey_dec_vq(bgds, tcds)[-1]
      covs.append(agg_time([vq_tot, vq_tcs, vq_t2b, vq_b2t, vq_bkg], latnm=ALAT, antisym=True, diff=False))

   lbls = ["$v'q'$", "$v_{TC}'q_{TC}'$", "$v_{TC}'q_{bg}'$", "$v_{bg}'q_{TC}'$", "$v_{bg}'q_{bg}'$"]

   #eddy_terms = agg_time([ls[-1] for ls in [vq_tot, vq_tcs, vq_bkg, vq_t2b, vq_b2t]], antisym=True, diff=False)
   
   #[plt.plot(et.latitudes, et, label=lbls[ii]) for ii, et in enumerate(eddy_terms)]
   #plt.plot(eddy_terms[1].latitudes, eddy_terms[1] + eddy_terms[3] + eddy_terms[4], c='black')
   #plt.legend()
   #plt.savefig('vpqp_covar_terms_test.png')
   #plt.show()

   plt.rc('font', size=16)
   fig, axes = plt.subplots(2, 3, figsize=(22, 9), sharex=True, sharey=True)

   for ii, ax in enumerate(axes.ravel()):
       if ii == 4: continue
       ixh = IXHORS[ii]
       ax.plot(YSCL(covs[ixh][0][ALAT]), covs[ixh][0], label='mean')
       ax.plot(YSCL(covs[ixh][0][ALAT]), sum(covs[ixh][1:4]), label='TCs')

       ax.set_title(TTLS[ixh])
       ax.set_title('(%s)' % chr(ord('a') + ii), loc='left')
       ax.axhline(0, c='gray')
       [ax.axvline(yl, c='gray', lw=0.5) for yl in YLOC]
       ax.set_xticks(YLOC, YLAB)
       ax.tick_params(top=True, right=True, labelbottom=True, labelleft=True)
       ax.set_xlabel('Latitude [°]')
       ax.set_ylabel('$\\tau_x$' + ('' if ii == 1 else ' anomaly') + ' [N m$^{-2}$]')
       ax.legend(loc='upper left')

       #if not ixh == CTLIX:
       #   ax.set_xlim(YSCL(tcsplt[ixh][ALAT][0]), YSCL(tcsplt[ixh][ALAT][-1]))
       #   ax.set_ylim(-.015, .015)

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
   plt.savefig(f'TC_masked_VQ850.pdf', bbox_inches='tight')
   plt.show()
   plt.close()

if __name__ == '__main__':
   main()
