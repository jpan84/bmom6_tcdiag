import numpy as np
import xarray as xr
from sznl_funcs import stack_hemi_sznl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gspec

TCDENS = '/glade/u/home/jpan/aquaptc/tempest/250725_density_sznl/tcdens.nc' #output from track_dens.py
OCN_YZ = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/ocn_yz_plts_diff/sznlzm.nc' #output from plt_ocn_yz.py
H1IALL = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/tcfields2mps_%s_full/means_all.nc' #output from tcfieldcontrib1.py
H1ITCS = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/tcfields2mps_%s_full/means_tcs.nc'
HU_HTF = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/sznl_transports_newdiff/HT_HU.nc' #output from sznl_transports_1.py

ALIS = ['250415_unseed', '250417_ctrl', '250416_seed1x1']

ZSCL = np.vectorize(lambda zl: zl / 2850 if zl <= 2000 else 2000 / 2850 + (zl - 2000) / 2000 * 850 / 2850)
ZLAB = np.arange(1000, 5000, 1000)
ZLOC = ZSCL(ZLAB)

LIMLAB = np.arange(200, 1500, 200)
LIMLOC = ZSCL(LIMLAB)

YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)

def main():
   tcdens = xr.open_dataset(TCDENS)
   ocn_yz = xr.open_dataset(OCN_YZ)
   h1iall = xr.concat([xr.open_dataset(H1IALL % al).expand_dims(case=[al]) for al in ALIS], dim='case')
   h1itcs = xr.concat([xr.open_dataset(H1ITCS % al).expand_dims(case=[al]) for al in ALIS], dim='case')
   huhtds = xr.open_dataset(HU_HTF)

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
         tauxall = fldnet(h1iall, 'TAUX').isel(case=jj) if jj == 1 else fldnet(h1iall, 'TAUX').isel(case=jj) - fldnet(h1iall, 'TAUX').isel(case=1)
         tauxall = tauxall.sel(season=szn)
         #tauxtcs = fldnet(h1itcs, 'TAUX').sel(season=szn).isel(case=jj)

         ax_top.plot(YSCL(tcdens['lat']), pltdens, linestyle='dashed', color='black')
         acelim = (-5, 5) if jj <= 1 else (-25, 25)
         ax_top.set_ylim(*acelim)
         ax_top.set_ylabel('ACE')
         axt2 = ax_top.twinx()
         axt2.plot(YSCL(h1iall['latitudes']), tauxall, color='purple')
         taulim = (-.25, .25) if jj == 1 else (-.05, .05)
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
   sf_ctkwargs = {'colors': 'green', 'levels': 2.**np.arange(-2, 7, 1)}
   sf_ctkwargs['levels'] = np.concatenate((-sf_ctkwargs['levels'][::-1], sf_ctkwargs['levels']))
   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (24, 9)
   fig, axes = plt.subplots(2, 3, layout='constrained')
   fig.suptitle('Ocean Heat Uptake [W m$^{-2}$], ACE [10$^4$ kt$^2$ yr$^{-1}$ (10$^6$ km$^2$)$^{-1}$]\
                  \nOcean Eulerian Mean Mass Streamfunction, Potential Temp [K]')
   for ii, szn in enumerate(ocn_yz['season']):
      for jj in range(len(ALIS)):
         outer_spec = axes[ii, jj].get_subplotspec()
         gs = gspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_spec, height_ratios=[1, 4])
         ax_top = fig.add_subplot(gs[0])
         ax_bot = fig.add_subplot(gs[1], sharex=ax_top)
         axes[ii, jj].axis('off')

         sf = (ocn_yz['vmo']).isel(case=jj)
         if jj != 1:
            sf = sf - (ocn_yz['vmo']).isel(case=1)
         sf = sf.sel(season=szn)
         expo = 1e10 if jj == 1 else 1e9
         ax_bot.contour(YSCL(ocn_yz['yq']), ZSCL(ocn_yz['zl']), sf / expo, **sf_ctkwargs)
         ax_bot.set_xticks(YLOC, YLAB)
         ax_bot.set_xlim(-.75, .75)
         ax_bot.set_yticks(LIMLOC, LIMLAB)
         ax_bot.set_ylim(ZSCL(800), ZSCL(0))
         ax_bot.set_xlabel('Latitude [°]')
         ax_bot.set_ylabel('Depth [m]')

         thta = ocn_yz['thetao'].isel(case=jj) if jj == 1 else ocn_yz['thetao'].isel(case=jj) - ocn_yz['thetao'].isel(case=1)
         th_kw = dict(cmap='cividis', levels=np.arange(0, 37, 4))
         if jj != 1:
            th_kw = dict(cmap='RdBu_r', levels=np.arange(-.5, .51, .05))
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
         ax_top.set_ylim(*acelim)
         ax_top.set_ylabel('ACE')
         axt2 = ax_top.twinx()
         axt2.plot(YSCL(huhtds['lat']), ohuall, color='maroon')
         ohulim = (-150, 150) if jj == 1 else (-20, 20)
         axt2.set_ylim(*ohulim)
         axt2.set_ylabel('OHU', color='maroon')
         axt2.tick_params(axis='y', colors='maroon')
         #axt2.plot(YSCL(h1itcs['latitudes']), tauxtcs, linestyle='-.')
         axt2.hlines(0, -1, 1, linestyles='dotted', colors='gray')

   #fig.tight_layout()
   fig.get_layout_engine().set(w_pad=.25, h_pad=0.15)
   #plt.show()
   plt.savefig('ocn_yz_plts_diff/ACE_OHU_vmo_theto.png', bbox_inches='tight')
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
