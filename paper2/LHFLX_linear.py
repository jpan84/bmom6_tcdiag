import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post')
from paths import ARCHRT, ALIA, CTLIX, CASENAMES
import consts as c
import xarray as xr
from sznl_funcs import stack_hemi_sznl, monthly2sznl, agg_time
import matplotlib.pyplot as plt

zm2d = '/glade/campaign/univ/upsu0032/jpan_aquaptc/%s/atm/uxzm_hist_h0a_noncons_-60.0_60.0_1.5_LHFLX_U10_TREFHT_QREFHT_TS_PS_TAUX_TAUY.nc'
zm3d = '/glade/campaign/univ/upsu0032/jpan_aquaptc/%s/atm/uxzm_hist_h0a_noncons_-60.0_60.0_1.5_U_V_UU_VV.nc'

def main():
   dss2d = xr.concat([xr.open_dataset(zm2d % cs).expand_dims(case=[ALIA[ii]]) for ii, cs in enumerate(CASENAMES)], dim='case')
   dss2d = dss2d.assign(MAGTAU=(dss2d['TAUX']**2 + dss2d['TAUY']**2)**0.5)
   dss3d = xr.concat([xr.open_dataset(zm3d % cs).expand_dims(case=[ALIA[ii]]) for ii, cs in enumerate(CASENAMES)], dim='case')
   dss3d = dss3d.assign(WSPD=(dss3d['UU'] + dss3d['VV'])**0.5)
   dss3d = dss3d.assign(KE=0.5 * (dss3d['UU'] + dss3d['VV']))
   dss3d = dss3d.assign(MKE=0.5 * (dss3d['U']**2 + dss3d['V']**2))
   dss3d = dss3d.assign(EKE=dss3d['KE'] - dss3d['MKE'])

   print(dss2d)
   print(dss3d)

   hy2d = dss2d.map(agg_time)
   hy3d = dss3d.map(agg_time).drop_vars('V')
   hy_v = agg_time(dss3d['V'], antisym=True)

   print(hy2d)

   dif2d = (hy2d - hy2d.isel(case=CTLIX)).drop_sel(case=ALIA[CTLIX])
   dif3d = (hy3d - hy3d.isel(case=CTLIX)).drop_sel(case=ALIA[CTLIX])

   frac2d = dif2d / hy2d.isel(case=CTLIX)
   frac3d = dif3d / hy3d.isel(case=CTLIX)

   frac_mke = dif3d['MKE'] / hy3d['KE'].isel(case=CTLIX)
   frac_eke = dif3d['EKE'] / hy3d['KE'].isel(case=CTLIX)

   #ke = 0.5 * (hy3d['UU'] + hy3d['VV']).isel(lev=-1)
   #mke = 0.5 * (hy3d['U']**2 + hy_v**2).isel(lev=-1)
   #eke = ke - mke

   #ke_frac = (ke - ke.isel(case=CTLIX)).drop_isel(case=CTLIX) / ke.isel(case=CTLIX)
   #mke_frac = (mke - mke.isel(case=CTLIX)).drop_isel(case=CTLIX) / mke.isel(case=CTLIX)
   #eke_frac = (eke - eke.isel(case=CTLIX)).drop_isel(case=CTLIX) / eke.isel(case=CTLIX)

   plt.rcParams['figure.figsize'] = (16, 4)
   fig, axes = plt.subplots(1, 4)

   for ii, ax in enumerate(axes):
      ax.plot(dif2d.latitudes, frac2d['LHFLX'].isel(case=ii), lw=3, c='black', label='LHF')
      ax.plot(dif2d.latitudes, frac2d['U10'].isel(case=ii), lw=1, c='green', label='U10')
      #ax.plot(dif2d.latitudes, frac3d['WSPD'].isel(case=ii, lev=-1), lw=2, c='blue', label='UBOT')
      ax.plot(dif2d.latitudes, 0.5 * frac2d['MAGTAU'].isel(case=ii), lw=1, c='red', label='TAU')
      #ax.plot(dif3d.latitudes, 0.5 * (frac_eke.isel(case=ii, lev=-1) + frac_mke.isel(case=ii, lev=-1)), c='purple', label='total KE')
      #ax.plot(dif3d.latitudes, 0.5 * frac_mke.isel(case=ii, lev=-1), c='purple', ls='dashed', label='MKE')
      #ax.plot(dif3d.latitudes, 0.5 * frac_eke.isel(case=ii, lev=-1), c='purple', ls='dotted', label='EKE')
      plt.legend()
      ax.axhline(0, lw=0.5, c='gray')
      #ax.set_xlim(5, 40)

   fig.tight_layout()
   plt.show()

if __name__ == '__main__':
   main()
