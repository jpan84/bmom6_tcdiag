import numpy as np
import xarray as xr
from sznl_funcs import stack_hemi_sznl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gspec
import matplotlib.ticker as mticker

TCDENS = '/glade/u/home/jpan/aquaptc/tempest/250725_density_sznl/tcdens.nc' #output from track_dens.py
OCN_YZ = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/ocn_yz_plts_diff/sznlzm.nc' #output from plt_ocn_yz.py
H1IALL = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/tcfields2mps_%s_full/means_all.nc' #output from tcfieldcontrib1.py
H1ITCS = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/tcfields2mps_%s_full/means_tcs.nc'
HU_HTF = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/sznl_transports_newdiff/HT_HU.nc' #output from sznl_transports_1.py
CAM_2D = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/linevslat_new_sznl_diff/sznlzm.nc' #output from sznl_zm_ux_1.py
ATMTEM = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/streamf_sznl_0005-0015/%s'

TEMNCS = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl_b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed_TEM.nc',\
          'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl__TEM.nc',
          'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl_b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1_TEM.nc']

ALIS = ['250415_unseed', '250417_ctrl', '250416_seed1x1']

ZSCL = np.vectorize(lambda zl: zl / 2850 if zl <= 2000 else 2000 / 2850 + (zl - 2000) / 2000 * 850 / 2850)
ZLAB = np.arange(1000, 5000, 1000)
ZLOC = ZSCL(ZLAB)

LIMLAB = np.arange(200, 1500, 200)
LIMLOC = ZSCL(LIMLAB)

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.arange(-60, 61, 10) #np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)

def main():
   tcdens = xr.open_dataset(TCDENS)
   ocn_yz = xr.open_dataset(OCN_YZ)
   h1iall = xr.concat([xr.open_dataset(H1IALL % al).expand_dims(case=[al]) for al in ALIS], dim='case')
   h1itcs = xr.concat([xr.open_dataset(H1ITCS % al).expand_dims(case=[al]) for al in ALIS], dim='case')
   huhtds = xr.open_dataset(HU_HTF)
   cam_2d = xr.open_dataset(CAM_2D)
   temds = xr.concat([xr.open_dataset(ATMTEM % TEMNCS[ii]).expand_dims(case=[al]) for ii, al in enumerate(ALIS)], dim='case')

   #mirror and stack the seasons of TC density
   for dv in tcdens.data_vars:
      stacked = stack_hemi_sznl(tcdens[dv], antisym=False, latnm='lat')
      tcdens = tcdens.assign(variables={dv: stacked})
   tcdens = tcdens.sel(season=['JJA', 'SON'])

   print(tcdens)
   print(tcdens['lat'])
   print(ocn_yz)
   print(h1iall)
   print(h1itcs)
   print(huhtds)

   ############### Plot TAUX and ACE density above vmo ####################################
   sf_cfkwargs = {'levels': 2.**np.arange(-2, 7, 1), 'norm': colors.SymLogNorm(2.**-1)}
   sf_cfkwargs['levels'] = np.concatenate((-sf_cfkwargs['levels'][::-1], sf_cfkwargs['levels']))
   sf_ctkwargs = {'colors': 'black', 'levels': 2.**np.arange(-2, 7, 1)}
   sf_ctkwargs['levels'] = np.concatenate((-sf_ctkwargs['levels'][::-1], sf_ctkwargs['levels']))
   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (24, 9)
   fig, axes = plt.subplots(2, 3, layout='constrained')
   fig.suptitle('Surface zonal stress [Pa m$^{-2}$], ACE [10$^4$ kt$^2$ yr$^{-1}$ (10$^6$ km$^2$)$^{-1}$]\
                  \nOcean Eulerian Mean Mass Streamfunction')
   for ii, szn in enumerate(ocn_yz['season']):
      for jj in range(len(ALIS)):
         outer_spec = axes[ii, jj].get_subplotspec()
         gs = gspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_spec, height_ratios=[1, 4])
         ax_top = fig.add_subplot(gs[0])
         ax_bot = fig.add_subplot(gs[1], sharex=ax_top)
         axes[ii, jj].axis('off')

         sf = ocn_yz['vmo'].isel(case=jj) if jj == 1 else ocn_yz['vmo'].isel(case=jj) - ocn_yz['vmo'].isel(case=1)
         sf = sf.sel(season=szn)
         expo = 1e10 if jj == 1 else 1e9
         csf = ax_bot.contourf(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / expo, cmap='coolwarm' if jj == 1 else 'bwr', **sf_cfkwargs)
         cax = ax_bot.inset_axes([1, 0, .05, 1], transform=ax_bot.transAxes)
         plt.colorbar(csf, label='[10$^{%d}$ kg/s]' % int(np.log10(expo)), cax=cax, pad=.01)
         ax_bot.contour(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / expo, **sf_ctkwargs)
         ax_bot.set_xticks(YLOC, YLAB)
         ax_bot.set_xlim(-1, 1)
         ax_bot.set_yticks(ZLOC, ZLAB)
         ax_bot.invert_yaxis()
         ax_bot.set_xlabel('Latitude [°]')
         ax_bot.set_ylabel('Depth [m]')

         #TODO: make sure TAUX includes all 10 years
         fldnet = lambda ds, fld: ds['%s_neg' % fld] + ds['%s_pos' % fld]
         pltdens = tcdens['ace'].isel(run=jj) if jj == 1 else tcdens['ace'].isel(run=jj) - tcdens['ace'].isel(run=1)
         pltdens = pltdens.sel(season=szn)
         tauxall = -cam_2d['TAUX'].isel(case=jj)# if jj == 1 else -cam_2d['TAUX'].isel(case=jj) + cam_2d['TAUX'].isel(case=1)
         #tauxall = fldnet(h1iall, 'TAUX').isel(case=jj) if jj == 1 else fldnet(h1iall, 'TAUX').isel(case=jj) - fldnet(h1iall, 'TAUX').isel(case=1)
         tauxall = tauxall.sel(season=szn)
         #tauxtcs = fldnet(h1itcs, 'TAUX').sel(season=szn).isel(case=jj)

         ax_top.plot(YSCL(tcdens['lat']), pltdens, linestyle='dashed', color='black')
         acelim = (-5, 5) if jj <= 1 else (-25, 25)
         ax_top.set_ylim(*acelim)
         ax_top.set_ylabel('ACE')
         axt2 = ax_top.twinx()
         axt2.plot(YSCL(cam_2d['latitudes']), tauxall, color='purple')
         taulim = (-.25, .25) if jj == 1 else (-.02, .02)
         axt2.set_ylim(*taulim)
         axt2.set_ylabel('$\\tau_x$', color='purple')
         axt2.tick_params(axis='y', colors='purple')
         #axt2.plot(YSCL(h1itcs['latitudes']), tauxtcs, linestyle='-.')
         axt2.hlines(0, -1, 1, linestyles='dotted', colors='gray')

   #fig.tight_layout()
   fig.get_layout_engine().set(w_pad=.25, h_pad=0.15)
   #plt.show()
   plt.savefig('ocn_yz_plts_diff/ACE_TAUX_vmo.png', bbox_inches='tight')
   plt.close()
   #########################################################################################


   ############### Plot OHU and ACE density above vmo, thetao #################################
   sf_ctkwargs = {'colors': 'green', 'levels': 10.**np.arange(9, 12)}#np.array([0.5, 2, 8])}#2.**np.arange(-2, 7, 1)}
   sf_ctkwargs['levels'] = np.concatenate((-sf_ctkwargs['levels'][::-1], sf_ctkwargs['levels']))
   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (24, 9)
   fig, axes = plt.subplots(2, 3, layout='constrained')
   fig.suptitle('Ocean Heat Uptake [W m$^{-2}$], ACE [10$^4$ kt$^2$ yr$^{-1}$ (10$^6$ km$^2$)$^{-1}$]\
                  \nOcean Residual Mass Streamfunction, Potential Temp [K]')
   acelims = [(-5, 5), (-5, 5), (-30, 30)]
   ohulims = [(-5, 5), (-120, 120), (-25, 25)]
   thelims = [(-.2, .21, .04), (0, 37, 4), (-1.2, 1.3, .4)]

   for ii, szn in enumerate(ocn_yz['season']):
      for jj in range(len(ALIS)):
         outer_spec = axes[ii, jj].get_subplotspec()
         gs = gspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_spec, height_ratios=[1, 4])
         ax_top = fig.add_subplot(gs[0])
         ax_bot = fig.add_subplot(gs[1], sharex=ax_top)
         axes[ii, jj].axis('off')

         sf = (ocn_yz['vmo'] + ocn_yz['vhGM']).isel(case=jj)
         if jj != 1:
            sf = sf - (ocn_yz['vmo'] + ocn_yz['vhGM']).isel(case=1)
         sf = sf.sel(season=szn)
         expo = 1e10 if jj == 1 else 1e9
         ax_bot.contour(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf, **sf_ctkwargs)
         ax_bot.set_xticks(YLOC, YLAB)
         ax_bot.set_xlim(-1 / np.sqrt(2), 1 / np.sqrt(2))
         ax_bot.set_yticks(LIMLOC, LIMLAB)
         ax_bot.set_ylim(ZSCL(800), ZSCL(0))
         ax_bot.set_xlabel('Latitude [°]')
         ax_bot.set_ylabel('Depth [m]')

         thta = ocn_yz['thetao'].isel(case=jj) if jj == 1 else ocn_yz['thetao'].isel(case=jj) - ocn_yz['thetao'].isel(case=1)
         th_kw = dict(cmap='cividis', levels=np.arange(0, 37, 4))
         if jj != 1:
            th_kw = dict(cmap='RdBu_r', levels=np.arange(*thelims[jj]))
         csf = ax_bot.contourf(YSCL(ocn_yz['yh']), ZSCL(ocn_yz['zl']), thta.sel(season=szn), extend='both', **th_kw)
         cax = ax_bot.inset_axes([1, 0, .05, 1], transform=ax_bot.transAxes)
         plt.colorbar(csf, label='Potential temp [K]', cax=cax, pad=.01)

         pltdens = tcdens['ace'].isel(run=jj) if jj == 1 else tcdens['ace'].isel(run=jj) - tcdens['ace'].isel(run=1)
         pltdens = pltdens.sel(season=szn)
         ohuall = huhtds['ohu'].isel(case=jj) if jj == 1 else huhtds['ohu'].isel(case=jj) - huhtds['ohu'].isel(case=1)
         ohuall = ohuall.sel(season=szn)
         #tauxtcs = fldnet(h1itcs, 'TAUX').sel(season=szn).isel(case=jj)

         ax_top.plot(YSCL(tcdens['lat']), pltdens, linestyle='dashed', color='black')
         acelim = (-5, 5) if jj <= 1 else (-25, 25)
         ax_top.set_ylim(*acelims[jj])
         ax_top.set_ylabel('ACE')
         axt2 = ax_top.twinx()
         axt2.plot(YSCL(huhtds['lat']), ohuall, color='maroon')
         ohulim = (-120, 120) if jj == 1 else (-13, 13)
         axt2.set_ylim(*ohulims[jj])
         axt2.set_ylabel('OHU', color='maroon')
         axt2.tick_params(axis='y', colors='maroon')
         #axt2.plot(YSCL(h1itcs['latitudes']), tauxtcs, linestyle='-.')
         axt2.hlines(0, -1, 1, linestyles='dotted', colors='gray')

   #fig.tight_layout()
   fig.get_layout_engine().set(w_pad=.25, h_pad=0.15)
   #plt.show()
   plt.savefig('ocn_yz_plts_diff/ACE_OHU_resid_thetao.png', bbox_inches='tight')
   plt.close()
   #########################################################################################

   '''
   ############### Plot atmo Eulerian sf above vresid, thetao #################################
   sf_cfkwargs = {'levels': 2.**np.arange(-2, 7, 1), 'norm': colors.SymLogNorm(2.**-1)}
   sf_cfkwargs['levels'] = np.concatenate((-sf_cfkwargs['levels'][::-1], sf_cfkwargs['levels']))
   sf_ctkwargs = {'colors': 'black', 'levels': 2.**np.arange(-2, 7, 1)}
   sf_ctkwargs['levels'] = np.concatenate((-sf_ctkwargs['levels'][::-1], sf_ctkwargs['levels']))
   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (24, 12)
   # --- NEW/MODIFIED CODE START ---
   
   # Define the master GridSpec
   # 2 rows for seasons (ii), 6 columns for 3 cases (jj) + 3 colorbars
   # width_ratios: [1.0, 0.05] repeated 3 times.
   gs_ratios = [1.0, 0.05] * 3
   fig = plt.figure(figsize=plt.rcParams['figure.figsize'])
   gs = fig.add_gridspec(nrows=2, ncols=6, wspace=0.1, hspace=0.1,
                         height_ratios=[1, 1], # You only have 2 vertical plots per main cell, but we'll use subplots for separation.
                         width_ratios=gs_ratios)
   
   fig.suptitle('Atmosphere Eulerian mean mass streamfunction\
                               \nOcean Residual Mass Streamfunction')
   
   # Remove the old fig, axes = plt.subplots(...) line
   # fig, axes = plt.subplots(2, 3, layout='constrained') 
   
   # --- NEW/MODIFIED CODE END -
   #fig, axes = plt.subplots(2, 3, layout='constrained')
   #fig.suptitle('Atmosphere Eulerian mean mass streamfunction\
   #               \nOcean Residual Mass Streamfunction')
   for ii, szn in enumerate(ocn_yz['season']):
      for jj in range(len(ALIS)):
          
         outer_spec = axes[ii, jj].get_subplotspec()
         #gs = gspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_spec, height_ratios=[1, 1])
         gs = gspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer_spec, height_ratios=[1, 1], width_ratios=[1, 0.05], wspace=0.05) # Add a column for cax
         ax_top = fig.add_subplot(gs[0, 0]) # Change index to 0, 0
         ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top) # Change index to 1, 0
         cax = fig.add_subplot(gs[:, 1]) # Create the colorbar axis spanning both rows in the second column
         #ax_top = fig.add_subplot(gs[0])
         #ax_bot = fig.add_subplot(gs[1], sharex=ax_top)
         axes[ii, jj].axis('off')
         

         # --- MODIFIED CODE START ---            
         # Calculate the starting column index for this plot group (0, 2, or 4)
         col_start = jj * 2 
         
         # Use a sub-GridSpec *within the cell* to separate ax_top and ax_bot vertically
         # This is simpler and more reliable than trying to manage the vertical split 
         # with the top-level GridSpec.

         # Get the cell for the main plot area (spanning 2 rows, 1 column)
         plot_spec = gs[ii, col_start] 
         
         # Create a GridSpecFromSubplotSpec to manage the vertical split (ax_top/ax_bot)
         # This GridSpec now only manages the vertical split inside ONE of the main cells.
         gs_vertical = plot_spec.subgridspec(2, 1, height_ratios=[1, 1], hspace=0.0) 
         
         # Add the subplots
         ax_top = fig.add_subplot(gs_vertical[0, 0])
         ax_bot = fig.add_subplot(gs_vertical[1, 0], sharex=ax_top)
         
         # Add the colorbar axes to the next column (col_start + 1)
         cax = fig.add_subplot(gs[ii, col_start + 1]) 
         
         # The placeholder axes[ii, jj] and its .axis('off') are gone.
         
         # --- MODIFIED CODE END ---

         sf = (ocn_yz['vmo'] + ocn_yz['vhGM']).isel(case=jj)
         if jj != 1:
            sf = sf - (ocn_yz['vmo'] + ocn_yz['vhGM']).isel(case=1)
         sf = sf.sel(season=szn)
         expo = 1e10 if jj == 1 else 1e9
         csf = ax_bot.contourf(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / expo, cmap='coolwarm' if jj == 1 else 'bwr', **sf_cfkwargs)
         #cax = axes[ii][jj].inset_axes([1, 0, .05, 1], transform=axes[ii][jj].transAxes)
         plt.colorbar(csf, label='[10$^{%d}$ kg/s]' % int(np.log10(expo)), cax=cax, pad=.01)
         ax_bot.contour(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / expo, **sf_ctkwargs)
         ax_bot.set_xticks(YLOC, YLAB)
         ax_bot.set_xlim(-1, 1)
         ax_bot.set_yticks(ZLOC, ZLAB)
         ax_bot.invert_yaxis()
         ax_bot.set_xlabel('Latitude [°]')
         ax_bot.set_ylabel('Depth [m]')

         mmc = temds['PSI_EM'].isel(case=jj, season=ii)
         csf = ax_top.contourf(YSCL(temds['lat']), temds['plev'], mmc / expo, cmap='coolwarm' if jj == 1 else 'bwr', **sf_cfkwargs)
         ax_top.contour(YSCL(temds['lat']), temds['plev'], mmc / expo, **sf_ctkwargs)
         ax_top.set_ylim(1000, 100)
         ax_top.set_yscale('log')
         ax_top.yaxis.set_minor_formatter(mticker.ScalarFormatter())
         ax_top.yaxis.set_major_formatter(mticker.ScalarFormatter())
         ax_top.set_ylabel('Pressure [hPa]')
         ax_top.set_yticks(np.arange(100, 1001, 100))
         ax_top.set_xticklabels([])

   #fig.tight_layout()
   #fig.get_layout_engine().set(w_pad=.25, h_pad=0.15)
   #plt.show()
   plt.savefig('ocn_yz_plts_diff/atmo_ocn_sfunc.png', bbox_inches='tight')
   plt.close()
   #########################################################################################
   '''

   ############### Plot atmo Eulerian sf above vresid, thetao #################################
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
   # Increased gap ratio from 0.1 to 0.5 for wider separation between the 3 major columns.
   
   fig = plt.figure(figsize=plt.rcParams['figure.figsize'])
   gs = fig.add_gridspec(nrows=2, ncols=8, wspace=0.1, hspace=0.2,
                         height_ratios=[1, 1],
                         width_ratios=gs_ratios)
   
   fig.suptitle('Atmosphere Eulerian mean mass streamfunction\
                   \nOcean Residual Mass Streamfunction')
   
   # --- NEW/MODIFIED CODE END ---
   
   for ii, szn in enumerate(ocn_yz['season']):
      for jj in range(len(ALIS)):
         
         # Calculate the starting column index for this plot group (0, 3, or 6)
         # 0: First plot group (Cols 0, 1)
         # 1: Second plot group (Cols 3, 4)
         # 2: Third plot group (Cols 6, 7)
         col_start = jj * 3 

         # 1. Create sub-GridSpec for vertical plot splitting (ax_top/ax_bot)
         # The main plot area is in the first column of the group (col_start)
         plot_spec = gs[ii, col_start] 
         
         # Create a GridSpecFromSubplotSpec to manage the tight vertical split (ax_top/ax_bot)
         gs_vertical = plot_spec.subgridspec(2, 1, height_ratios=[1, 1], hspace=0.0)
         
         # Add the subplots
         ax_top = fig.add_subplot(gs_vertical[0, 0])
         ax_bot = fig.add_subplot(gs_vertical[1, 0], sharex=ax_top)
         
         # Add the colorbar axes to the second column of the group (col_start + 1)
         cax = fig.add_subplot(gs[ii, col_start + 1]) 
         
         # --- NEW: Set Thicker Borders (Spines) for ax_top and ax_bot ---
         border_thickness = 2.0
         for axis in [ax_top, ax_bot]:
            for spine in axis.spines.values():
                spine.set_linewidth(border_thickness)
         # ---------------------------------------------------------------
         
         sf = (ocn_yz['vmo'] + ocn_yz['vhGM']).isel(case=jj)
         if jj != 1:
            sf = sf - (ocn_yz['vmo'] + ocn_yz['vhGM']).isel(case=1)
         sf = sf.sel(season=szn)
         expo = 1e10 if jj == 1 else 1e9
         
         # Determine correct symmetric ticks for the current expo
         current_ticks = base_sym_ticks# * (1e10 / expo)

         csf = ax_bot.contourf(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / expo, cmap='coolwarm' if jj == 1 else 'bwr', **sf_cfkwargs)
         
         # --- MODIFIED: Use the newly created cax and defined ticks ---
         plt.colorbar(csf, label='[10$^{%d}$ kg/s]' % int(np.log10(expo)), 
                      cax=cax, pad=.01, ticks=current_ticks)
         # -------------------------------------------------------------
         
         ax_bot.contour(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / expo, **sf_ctkwargs)
         ax_bot.set_xticks(YLOC, YLAB) # Keeps xtick labels on ax_bot
         ax_bot.set_xlim(-1, 1)
         ax_bot.set_yticks(ZLOC, ZLAB)
         ax_bot.invert_yaxis()
         ax_bot.set_xlabel('Latitude [°]')
         ax_bot.set_ylabel('Depth [m]')

         mmc = temds['PSI_EM'].isel(case=jj, season=ii)
         csf = ax_top.contourf(YSCL(temds['lat']), temds['plev'], mmc / expo, cmap='coolwarm' if jj == 1 else 'bwr', extend='both', **sf_cfkwargs)
         ax_top.contour(YSCL(temds['lat']), temds['plev'], mmc / expo, **sf_ctkwargs)
         ax_top.set_ylim(1000, 100)
         ax_top.set_yscale('log')
         ax_top.yaxis.set_minor_formatter(mticker.ScalarFormatter())
         ax_top.yaxis.set_major_formatter(mticker.ScalarFormatter())
         ax_top.set_ylabel('Pressure [hPa]')
         ax_top.set_yticks(np.arange(100, 1001, 100))

         label_index = (ii * 3) + jj
         panel_label = f'({chr(97 + label_index)})'
         ax_top.set_title(panel_label, loc='left')
         
         # --- NEW/MODIFIED: Ensure no xtick labels on ax_top ---
         ax_top.tick_params(labelbottom=False) 
         # ax_top.set_xticklabels([]) is now redundant but keeping tick_params is cleaner
         # -------------------------------------------------------

   # 231    #fig.tight_layout() # tight_layout is not used with GridSpec
   # 232    #fig.get_layout_engine().set(w_pad=.25, h_pad=0.15) # Remove constrained_layout engine commands
   # 233    #plt.show()
   
   # Assuming the output commands remain the same
   plt.savefig('ocn_yz_plts_diff/atmo_ocn_sfunc.png', bbox_inches='tight')
   plt.close()
   #########################################################################################

   ############### Plot OHT and ACE density above vmo_resid, thetao #################################
   sf_ctkwargs = {'colors': 'green', 'levels': 2.**np.arange(-2, 7, 1)}
   sf_ctkwargs['levels'] = np.concatenate((-sf_ctkwargs['levels'][::-1], sf_ctkwargs['levels']))
   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (24, 9)
   fig, axes = plt.subplots(2, 3, layout='constrained')
   fig.suptitle('Ocean Heat Transport [W], ACE [10$^4$ kt$^2$ yr$^{-1}$ (10$^6$ km$^2$)$^{-1}$]\
                  \nOcean Residual Mass Streamfunction (not difference), Potential Temp [K]')
   for ii, szn in enumerate(ocn_yz['season']):
      for jj in range(len(ALIS)):
         outer_spec = axes[ii, jj].get_subplotspec()
         gs = gspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_spec, height_ratios=[1, 4])
         ax_top = fig.add_subplot(gs[0])
         ax_bot = fig.add_subplot(gs[1], sharex=ax_top)
         axes[ii, jj].axis('off')

         sf = (ocn_yz['vmo'] + ocn_yz['vhGM']).isel(case=jj)
         #if jj != 1:
         #   sf = sf - (ocn_yz['vmo']).isel(case=1)
         sf = sf.sel(season=szn)
         expo = 1e10 #if jj == 1 else 1e9
         ax_bot.contour(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / expo, **sf_ctkwargs)
         ax_bot.set_xticks(YLOC, YLAB)
         ax_bot.set_xlim(-1, 1)
         ax_bot.set_yticks(ZLOC, ZLAB)
         ax_bot.invert_yaxis()
         ax_bot.set_xlabel('Latitude [°]')
         ax_bot.set_ylabel('Depth [m]')

         thta = ocn_yz['thetao'].isel(case=jj) if jj == 1 else ocn_yz['thetao'].isel(case=jj) - ocn_yz['thetao'].isel(case=1)
         th_kw = dict(cmap='cividis', levels=np.arange(0, 37, 4))
         if jj != 1:
            th_kw = dict(cmap='RdBu_r', levels=np.arange(-.2, .21, .02))
         csf = ax_bot.contourf(YSCL(ocn_yz['yh']), ZSCL(ocn_yz['zl']), thta.sel(season=szn), extend='both', **th_kw)
         cax = ax_bot.inset_axes([1, 0, .05, 1], transform=ax_bot.transAxes)
         plt.colorbar(csf, label='Potential temp [K]', cax=cax, pad=.01)

         pltdens = tcdens['ace'].isel(run=jj) if jj == 1 else tcdens['ace'].isel(run=jj) - tcdens['ace'].isel(run=1)
         pltdens = pltdens.sel(season=szn)
         ohtall = huhtds['oht'].isel(case=jj) if jj == 1 else huhtds['oht'].isel(case=jj) - huhtds['oht'].isel(case=1)
         ohtall = ohtall.sel(season=szn)
         #tauxtcs = fldnet(h1itcs, 'TAUX').sel(season=szn).isel(case=jj)

         ax_top.plot(YSCL(tcdens['lat']), pltdens, linestyle='dashed', color='black')
         acelim = (-5, 5) if jj <= 1 else (-25, 25)
         ax_top.set_ylim(*acelim)
         ax_top.set_ylabel('ACE')
         axt2 = ax_top.twinx()
         axt2.plot(YSCL(huhtds['yq']), ohtall, color='blue')
         ohtlim = (-8e15, 8e15) if jj == 1 else (-4e14, 4e14)
         axt2.set_ylim(*ohtlim)
         axt2.set_ylabel('OHT', color='blue')
         axt2.tick_params(axis='y', colors='blue')
         #axt2.plot(YSCL(h1itcs['latitudes']), tauxtcs, linestyle='-.')
         axt2.hlines(0, -1, 1, linestyles='dotted', colors='gray')

   #fig.tight_layout()
   fig.get_layout_engine().set(w_pad=.25, h_pad=0.15)
   #plt.show()
   plt.savefig('ocn_yz_plts_diff/ACE_OHT_resid_theto.png', bbox_inches='tight')
   plt.close()
   #########################################################################################

if __name__ == '__main__':
   main()
