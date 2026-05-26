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

   ############### Plot atmo MSE,DSE above ocean theta,MLdepth #################################
   #TODO: set contour levels
   sf_cfkwargs = {'levels': 2.**np.arange(-2, 7, 1), 'norm': colors.SymLogNorm(2.**-1)}
   sf_cfkwargs['levels'] = np.concatenate((-sf_cfkwargs['levels'][::-1], sf_cfkwargs['levels']))
   sf_ctkwargs = {'colors': 'black', 'levels': 2.**np.arange(-2, 7, 1)}
   sf_ctkwargs['levels'] = np.concatenate((-sf_ctkwargs['levels'][::-1], sf_ctkwargs['levels']))
   
   # Define symmetric ticks for the colorbars
   pos_ticks = 2.**np.arange(0, 7, 2) # [1., 4., 16., 64.]
   base_sym_ticks = np.concatenate((-pos_ticks[::-1], [0], pos_ticks))
   
   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (24, 12)
   
   # --- NEW/MODIFIED CODE START: Unified GridSpec for spacing control ---
   
   # Total Grid Dimensions: 2 rows, 8 columns (3 Plot, 3 Cax, 2 Gap)
   # [Plot, Cax, Gap, Plot, Cax, Gap, Plot, Cax]
   gs_ratios = [1.3, 0.05, 0.5, 1.3, 0.05, 0.5, 1.3, 0.05] 
   
   fig = plt.figure(figsize=plt.rcParams['figure.figsize'])
   gs = fig.add_gridspec(nrows=2, ncols=8, wspace=0.1, hspace=0.2,
                         height_ratios=[1, 1],
                         width_ratios=gs_ratios)
   
   #fig.suptitle('Atmosphere Eulerian mean mass streamfunction\
                   #\nOcean Residual Mass Streamfunction')
   
   # --- NEW/MODIFIED CODE END ---
   
   #for ii, szn in enumerate(ocn_yz['season']):
   #   for jj in range(len(ALIS)):

   for rr in range(gs.nrows):
      for cc in range(3): #remove the hardcoded 3
         ix1d = rr * 3 + cc
         ixh = IXHORS[ix1d]
         col_start = cc * 3
         sp = gs[rr, col_start]

         # Create a GridSpecFromSubplotSpec to manage the tight vertical split (ax_top/ax_bot)
         gs_vertical = sp.subgridspec(2, 1, height_ratios=[1, 0.5], hspace=0.0)
         
         # Add the subplots
         ax_top = fig.add_subplot(gs_vertical[0, 0])
         ax_bot = fig.add_subplot(gs_vertical[1, 0], sharex=ax_top)
         
         # Add the colorbar axes to the second column of the group (col_start + 1)
         cax = fig.add_subplot(gs[rr, col_start + 1]) 
         
         # --- NEW: Set Thicker Borders (Spines) for ax_top and ax_bot ---
         border_thickness = 2.0
         for axis in [ax_top, ax_bot]:
            for spine in axis.spines.values():
                spine.set_linewidth(border_thickness)
         # ---------------------------------------------------------------
         
   
         # Determine correct symmetric ticks for the current expo
         current_ticks = base_sym_ticks# * (1e10 / expo)
   
         csf = ax_bot.contourf(YSCL(othplt[ixh]['lat']), othplt[ixh]['zl'], othplt[ixh], cmap='magma' if ixh == CTLIX else 'bwr')#, **sf_cfkwargs)
         
         # --- MODIFIED: Use the newly created cax and defined ticks ---
         #plt.colorbar(csf, label='[10$^{%d}$ kg/s]' % int(np.log10(expo)), 
         #             cax=cax, pad=.01, ticks=current_ticks)
         # -------------------------------------------------------------
         
         #ax_bot.contour(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / expo, **sf_ctkwargs)
         ax_bot.plot(YSCL(omlplt[ixh]['lat']), omlplt[ixh])
         ax_bot.set_xticks(YLOC, ['' if yl % 20 else yl for yl in YLAB]) # Keeps xtick labels on ax_bot
         ax_bot.set_xlim(-1, 1)
         ax_bot.set_ylim(0, 600)
         #ax_bot.set_yticks(ZLOC, ZLAB)
         ax_bot.invert_yaxis()
         ax_bot.set_xlabel('Latitude [°]')
         ax_bot.set_ylabel('Depth [m]')
         ax_bot.tick_params(right=True)
         ax_bot.tick_params(length=7, width=1.3, which='both')

         csf = ax_top.contourf(YSCL(mseplt[ixh]['latitudes']), mseplt[ixh]['lev'], mseplt[ixh])#, cmap='coolwarm' if jj == 1 else 'bwr', extend='both', **sf_cfkwargs)
         ax_top.contourf(YSCL(dseplt[ixh]['latitudes']), dseplt[ixh]['lev'], dseplt[ixh])
         #ax_top.contour(YSCL(temds['lat']), temds['plev'], temds['PSI_EM'].isel(case=1, season=ii) / cexpo, **sf_ctkwargs)
         ax_top.set_ylim(1000, 100)
         ax_top.set_yscale('log')
         ax_top.yaxis.set_minor_formatter(mticker.ScalarFormatter())
         ax_top.yaxis.set_major_formatter(mticker.ScalarFormatter())
         ax_top.set_ylabel('Pressure [hPa]')
         ax_top.set_yticks(np.arange(200, 1001, 200))
         ax_top.tick_params(top=True, right=True, length=7, width=1.3, which='both')
         ax_top.tick_params(which='minor', labelsize=0)
   
         label_index = (rr * 3) + cc
         panel_label = f'({chr(97 + label_index)})'
         ax_top.set_title(panel_label, loc='left')

         ax_top.set_title(TTLS[ixh], fontsize=24)
         
         # --- NEW/MODIFIED: Ensure no xtick labels on ax_top ---
         ax_top.tick_params(labelbottom=False) 
         # ax_top.set_xticklabels([]) is now redundant but keeping tick_params is cleaner
         # -------------------------------------------------------

   
   # Assuming the output commands remain the same
   plt.savefig('thermo_state_plt.svg', bbox_inches='tight')
   plt.close()
   #########################################################################################


if __name__ == '__main__':
   main()
