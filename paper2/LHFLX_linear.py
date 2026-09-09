import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post')
from paths import ARCHRT, ALIA, CTLIX, CASENAMES
import consts as c
import xarray as xr
from sznl_funcs import stack_hemi_sznl, monthly2sznl, agg_time
import matplotlib.pyplot as plt

zm2d = '/glade/campaign/univ/upsu0032/jpan_aquaptc/%s/atm/uxzm_hist_h0a_noncons_-60.0_60.0_1.5_LHFLX_U10_TREFHT_QREFHT.nc'
zm3d = '/glade/campaign/univ/upsu0032/jpan_aquaptc/%s/atm/uxzm_hist_h0a_noncons_-60.0_60.0_1.5_U_V_UU_VV.nc'

def main():
   dss2d = xr.concat([xr.open_dataset(zm2d % cs).expand_dims(case=[ALIA[ii]]) for ii, cs in enumerate(CASENAMES)], dim='case')
   dss3d = xr.concat([xr.open_dataset(zm3d % cs).expand_dims(case=[ALIA[ii]]) for ii, cs in enumerate(CASENAMES)], dim='case')

   print(dss2d)
   print(dss3d)

   hy2d = dss2d.map(agg_time)
   hy3d = dss3d.map(agg_time).drop_vars('V')

   print(hy2d)

   dif2d = (hy2d - hy2d.isel(case=CTLIX)).drop_sel(case=ALIA[CTLIX])
   dif3d = hy3d - hy3d.isel(case=CTLIX)

   frac2d = dif2d / hy2d.isel(case=CTLIX)

   plt.rcParams['figure.figsize'] = (13, 4)
   fig, axes = plt.subplots(1, 4)

   for ii, ax in enumerate(axes):
      ax.plot(dif2d.latitudes, frac2d['LHFLX'].isel(case=ii))
      ax.plot(dif2d.latitudes, frac2d['U10'].isel(case=ii))

   fig.tight_layout()
   plt.show()

if __name__ == '__main__':
   main()
