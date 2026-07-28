import os
import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/')
from sznl_funcs import stack_hemi_sznl, monthly2sznl
import xarray as xr
import matplotlib.pyplot as plt

DIRI = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm'
totfil = 'uxzm_hist_h1i_noncons_-90.0_90.0_2.0_TAUX_PRECT_V850_Q850_V850.Q850.nc'
tcsfil = 'uxzm_nff4mps_h1i_noncons_-90.0_90.0_2.0_TAUX_PRECT_V850_Q850_V850.Q850.nc'
bkgfil = 'uxzm_nff4mpsinvert_h1i_noncons_-90.0_90.0_2.0_TAUX_PRECT_V850_Q850_V850.Q850.nc'

def rey_dec_vq(vds, qds, vvar='V850', qvar='Q850'):
   vq = None
   if vds is qds:
      vq = vds[vvar + '_dot_' + qvar]
   else:
      vq = xr.zeros_like(vds[vvar])
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

def agg_time_1case(da, latnm='latitudes', antisym=False):
   ymonmean = da.groupby('time.month').mean()
   twoszns = stack_hemi_sznl(monthly2sznl(ymonmean), antisym=antisym, latnm=latnm)
   halfyr = twoszns.mean(dim='season')
   return halfyr

totds = xr.open_dataset(os.path.join(DIRI, totfil))
tcsds = xr.open_dataset(os.path.join(DIRI, tcsfil))
bkgds = xr.open_dataset(os.path.join(DIRI, bkgfil))

vq_tot = rey_dec_vq(totds, totds)
vq_tcs = rey_dec_vq(tcsds, tcsds)
vq_bkg = rey_dec_vq(bkgds, bkgds)
vq_t2b = rey_dec_vq(tcsds, bkgds)
vq_b2t = rey_dec_vq(bkgds, tcsds)

eddy_terms = agg_time([ls[-1] for ls in [vq_tot, vq_tcs, vq_bkg, vq_t2b, vq_b2t]], antisym=True, diff=False)

[plt.plot(et.latitudes, et) for et in eddy_terms]
plt.show()
