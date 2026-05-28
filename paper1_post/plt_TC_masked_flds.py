import numpy as np
import os
import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/')
from paths import ARCHRT, ALIA, CTLIX, IXHORS
import xarray as xr
from sznl_funcs import stack_hemi_sznl, monthly2sznl
import matplotlib.pyplot as plt

TOTFIL = 'atm/uxzm_h0a_noncons_-90_90_0.25_TAUX.nc'
TCSFIL = 'atm/uxzm_h1i_noncons_-90_90_0.25_TC_R4_TAUX_PRECT.nc'
ALAT = 'latitudes'

TTLS = [ALIA[ii] + '$–$CTL' if ii != CTLIX else 'CTL' for ii in range(len(ALIA))]

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-60, 61, 10) #np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)

totdss = [xr.decode_cf(xr.open_dataset(os.path.join(ar, TOTFIL))) for ar in ARCHRT]
tcsdss = [xr.decode_cf(xr.open_dataset(os.path.join(ar, TCSFIL))) for ar in ARCHRT]

def main(pltvar='TAUX'):
   print([ds.time for ds in totdss])
   def agg_time(listdas, latnm='latitudes', diff=True):
      ymonmean = [da.groupby('time.month').mean() for da in listdas]
      twoszns = [stack_hemi_sznl(monthly2sznl(ym), antisym=False, latnm=latnm) for ym in ymonmean]
      halfyr = [ts.mean(dim='season') for ts in twoszns]
      if diff:
         return [hy - halfyr[CTLIX] if ii != CTLIX else hy for ii, hy in enumerate(halfyr)]
      return halfyr

   totplt = agg_time([ds[pltvar] for ds in totdss], latnm=ALAT)
   tcsplt = agg_time([ds[pltvar] for ds in tcsdss], latnm=ALAT)


   plt.rc('font', size=16)
   fig, axes = plt.subplots(2, 3, figsize=(22, 9))

   for ii, ax in enumerate(axes.ravel()):
       ixh = IXHORS[ii]
       ax.plot(YSCL(totplt[ixh][ALAT]), totplt[ixh])
       ax.plot(YSCL(tcsplt[ixh][ALAT]), tcsplt[ixh])

   plt.show()
   #plt.savefig('thermo_state_plt.svg', bbox_inches='tight')
   plt.close()


if __name__ == '__main__':
   main()
