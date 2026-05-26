import numpy as np
import os
import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/')
from paths import ARCHRT, ALIA, CTLIX, IXHORS
import consts as c
import xarray as xr
from sznl_funcs import stack_hemi_sznl, monthly2sznl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gspec
import matplotlib.ticker as mticker

AFIL = 'atm/uxzm_cons_-90_90_0.25_Z3_T_Q.nc'
OFIL = 'ocn/hist/cdo_zm_thetao_oml_mlotst.nc'
ALAT, OLAT = 'latitudes', 'lat'

TTLS = [ALIA[ii] + '$–$CTL' if ii != CTLIX else 'CTL' for ii in range(len(ALIA))]

ZSCL = np.vectorize(lambda zl: zl / 2850 if zl <= 2000 else 2000 / 2850 + (zl - 2000) / 2000 * 850 / 2850)
ZLAB = np.arange(1000, 5000, 1000)
ZLOC = ZSCL(ZLAB)

LIMLAB = np.arange(200, 1500, 200)
LIMLOC = ZSCL(LIMLAB)

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-60, 61, 10) #np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)

def main():
   ads = [xr.open_dataset(os.path.join(ar, AFIL)) for ar in ARCHRT] #TODO: make sure uxarray zonal_mean() doesn't drop coords
   ods = [xr.open_dataset(os.path.join(ar, OFIL)) for ar in ARCHRT]
   ads = [ad.isel(time=slice(0, ods[ii]['time'].size)).assign_coords(time=ods[ii]['time']) for ii, ad in enumerate(ads)] #slice for now ATM and OCN time mismatch when still running
   ads = [ds.assign_coords(lev=xr.open_dataset('/glade/u/home/jpan/grids/L26_hyb.nc')['lev']) for ds in ads]

   dse = [c.g * ds['Z3'] + c.cp * ds['T'] for ds in ads]
   mse = [c.g * ds['Z3'] + c.cp * ds['T'] + c.lv * ds['Q'] for ds in ads]

   def agg_time(listdas, latnm='latitudes', diff=True):
      ymonmean = [da.groupby('time.month').mean() for da in listdas]
      twoszns = [stack_hemi_sznl(monthly2sznl(ym.squeeze()), antisym=False, latnm=latnm) for ym in ymonmean]
      halfyr = [ts.mean(dim='season') for ts in twoszns]
      if diff:
         return [hy - halfyr[CTLIX] if ii != CTLIX else hy for ii, hy in enumerate(halfyr)]
      return halfyr

   dseplt = agg_time(dse)
   mseplt = agg_time(mse)
   othplt = agg_time([ds['thetao'] for ds in ods], latnm=OLAT)
   omlplt = agg_time([ds['oml'] for ds in ods], latnm=OLAT, diff=False)

   mse_kw = dict(levels=np.arange(-2000, 2001, 200), cmap='seismic')
   dse_kw = dict(levels=mse_kw['levels'][mse_kw['levels'] != 0])
   ctl_kw = dict(levels=np.arange(3e5, 4.01e5, 5e3), cmap='viridis')
   olevs = np.arange(-.2, .21, .02)
   olevs_c = np.arange(0, 40, 4)

   plt.rc('font', size=16)
   fig = plt.figure(figsize=(22, 13)) # Marginally wider aspect ratio to accommodate shared side colorbars or clean baselines

   # --- NEW GRIDSPEC ARCHITECTURE ---
   # 2 Main Rows for data plots. 1 Tiny bottom row for shared anomaly colorbars.
   # 3 Columns for your experiments. 1 Tiny right column for the Control absolute colorbars.
   gs = fig.add_gridspec(nrows=3, ncols=4, 
                         height_ratios=[1, 1, 0.08], 
                         width_ratios=[1, 1, 1, 0.06],
                         wspace=0.25, hspace=0.3)

   # Placeholders to capture mappables for colorbars
   cf_atmo_ctl, cf_atmo_anom = None, None
   cf_ocn_ctl, cf_ocn_anom = None, None

   for rr in range(2): # Process 2 rows of subplots
      for cc in range(3): # Process 3 experiment columns
          ix1d = rr * 3 + cc
          ixh = IXHORS[ix1d]
          
          # Target the standard plot grid cell
          sp = gs[rr, cc]
          gs_vertical = sp.subgridspec(2, 1, height_ratios=[1, 0.5], hspace=0.05)

          ax_top = fig.add_subplot(gs_vertical[0, 0])
          ax_bot = fig.add_subplot(gs_vertical[1, 0], sharex=ax_top)

          # --- Handle Spines & Styling ---
          border_thickness = 2.0
          for axis in [ax_top, ax_bot]:
              for spine in axis.spines.values():
                  spine.set_linewidth(border_thickness)

          # --- Plot Ocean Data (Bottom Panel) ---
          is_control = (ixh == CTLIX)
          ocn_cmap = 'magma' if is_control else 'bwr'
          
          csf_ocn = ax_bot.contourf(YSCL(othplt[ixh]['lat']), othplt[ixh]['zl'], othplt[ixh], cmap=ocn_cmap, levels=olevs_c if is_control else olevs)
          
          # Cache mapping handles for later colorbar generation
          if is_control:
              cf_ocn_ctl = csf_ocn
          elif cf_ocn_anom is None:
              cf_ocn_anom = csf_ocn

          ax_bot.plot(YSCL(omlplt[ixh]['lat']), omlplt[ixh], color='black', linewidth=1.5)
          ax_bot.set_xticks(YLOC, ['' if yl % 20 else yl for yl in YLAB])
          ax_bot.set_xlim(-1, 1)
          ax_bot.set_ylim(0, 400)
          ax_bot.invert_yaxis()
          
          if rr == 1: # Only label the true bottom row X-axis
              ax_bot.set_xlabel('Latitude [°]')
          else:
              ax_bot.tick_params(labelbottom=False)
              
          ax_bot.set_ylabel('Depth [m]') if cc == 0 else ax_bot.tick_params(labelleft=False)

          # --- Plot Atmosphere Data (Top Panel) ---
          csf_atmo = ax_top.contourf(YSCL(mseplt[ixh]['latitudes']), mseplt[ixh]['lev'], mseplt[ixh], 
                                     extend='both', **(mse_kw if not is_control else ctl_kw))

          css_atmo = ax_top.contour(YSCL(dseplt[ixh]['latitudes']), dseplt[ixh]['lev'], dseplt[ixh], 
                                     colors='lime', levels=(dse_kw if not is_control else ctl_kw)['levels'])
          
          if is_control:
              cf_atmo_ctl = csf_atmo
          elif cf_atmo_anom is None:
              cf_atmo_anom = csf_atmo

          ax_top.set_ylim(1000, 100)
          ax_top.set_yscale('log')
          ax_top.yaxis.set_minor_formatter(mticker.ScalarFormatter())
          ax_top.yaxis.set_major_formatter(mticker.ScalarFormatter())
          
          ax_top.set_ylabel('Pressure [hPa]') if cc == 0 else ax_top.tick_params(labelleft=False)
          ax_top.set_yticks(np.arange(200, 1001, 200))
          ax_top.tick_params(labelbottom=False)

          # Annotations
          label_index = (rr * 3) + cc
          ax_top.set_title(f'({chr(97 + label_index)})', loc='left', fontweight='bold')
          ax_top.set_title(TTLS[ixh], fontsize=20, pad=10)

   # --- CLEAN COLORBAR ROUTING ---
   # 1. Absolute Control Colorbars (Stacked on the far right column)
   cax_atmo_ctl = fig.add_subplot(gs[0, 3])
   fig.colorbar(cf_atmo_ctl, cax=cax_atmo_ctl, label='CTL Atmo MSE [J/kg]')

   cax_ocn_ctl = fig.add_subplot(gs[1, 3])
   fig.colorbar(cf_ocn_ctl, cax=cax_ocn_ctl, label='CTL Ocean Temp [°C]')

   # 2. Shared Anomaly Colorbars (Spanning columns 0 to 2 along the bottom row)
   # Split the bottom slot into two horizontal tracks using subgridspec
   gs_bottom_cbars = gs[2, 0:3].subgridspec(1, 2, wspace=0.4)
   
   cax_atmo_anom = fig.add_subplot(gs_bottom_cbars[0, 0])
   fig.colorbar(cf_atmo_anom, cax=cax_atmo_anom, orientation='horizontal', label='Anom Atmo MSE Anomaly [J/kg]')

   cax_ocn_anom = fig.add_subplot(gs_bottom_cbars[0, 1])
   fig.colorbar(cf_ocn_anom, cax=cax_ocn_anom, orientation='horizontal', label='Anom Ocean Temp Anomaly [°C]')


   # Assuming the output commands remain the same
   plt.savefig('thermo_state_plt.svg', bbox_inches='tight')
   plt.close()
   #########################################################################################


if __name__ == '__main__':
   main()
