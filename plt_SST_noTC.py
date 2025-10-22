from sznl_funcs import stack_hemi_sznl
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

FILI = ['/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/cw_15d_sznlzm_UNSEED-CTRL.nc',
        '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/cw_15d_sznlzm.nc',
        '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/cw_15d_sznlzm_SEED-CTRL.nc']
DIRO = './SST_noTC'
TTL = ['UNSEED–CTRL', 'CTRL', 'SEED–CTRL']
FLDNM = 'tos'

LNCLR = ['blue', 'orange']
YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.array([-90, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)
ZLIM = [(-.18, .3), (-.18, .3), (-1.2, 2.0)]

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   dss = [xr.open_dataset(fl).squeeze(drop=True).groupby('time.season').mean()\
          for fl in FILI]
   [print(ds) for ds in dss]
   toplt = [stack_hemi_sznl(ds[FLDNM], latnm='lat') for ds in dss]
   origds = xr.open_dataset('./linevslat_new_sznl_diff/sznlzm.nc')

   plt.rc('font', size=15)
   plt.rcParams['figure.figsize'] = (18, 4)
   subplot_kw = dict(xlim=(-1, 1))
   fig, axes = plt.subplot_mosaic([['(a)', '(b)', '(c)']], layout='constrained', sharey=False, subplot_kw=subplot_kw)

   ax1 = None
   for ii, (lbl, ax) in enumerate(axes.items()):
      [ax.axvline(yl, linewidth=0.5, color='gray') for yl in YLOC]
      ax.set_xlabel('Latitude [°]')
      ax.set_ylabel('[°C]')
      for jj, szn in enumerate(toplt[ii].season):
         ax.plot(YSCL(toplt[ii]['lat']), toplt[ii].sel(season=szn), color=LNCLR[jj], linestyle='dashdot', label='%s cold wakes' % str(szn.data))
         if ii in {0, 2}:
            ax.plot(YSCL(origds['latitudes']), origds['TS'].isel(case=ii).sel(season=szn), color=LNCLR[jj], label='%s SST' % str(szn.data))
         else:
            ax1 = ax.twinx() if ax1 is None else ax1
            ax1.plot(YSCL(origds['latitudes']), origds['TS'].isel(case=ii).sel(season=szn) - 273.15, color=LNCLR[jj], label='%s SST' % str(szn.data))
            ax1.set_ylabel('CTRL SST [°C]')
         if ii == 0:
            ax.legend(loc='lower center', fontsize=11, ncol=2)
         ax.axhline(y=0, linewidth=0.5, color='gray')
         ax.set_xticks(YLOC, ['' if yl % 30 else yl for yl in YLAB])
         ax.set_ylim(*ZLIM[ii])
         ax.set_title(TTL[ii])
         ax.set_title(lbl, loc='left')

   plt.savefig(os.path.join(DIRO, 'delta_cw_SST.png'))
   plt.close()

if __name__ == '__main__':
   main()
