#Joshua Pan jvp5930@psu.edu May 2025
#Overlay TEM streamf, atmo ener, track density
#Use data on pressure levels and structured grid, interpolated using ncremap (hy2pres.sh)

#import paths as pt
DIRIO = './streamf/'

import sys
sys.path.append('/glade/u/home/jpan/aquaptc/aquapgrid/')
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import consts as c
import linevslat

#GLOBS = [dr[1] for dr in pt.HDIRS(0)]

print(sys.argv)
#DIRIDX = int(sys.argv[1])
DIFF = bool(int(sys.argv[1])) #False if len(sys.argv) < 3 else bool(sys.argv[2])

SFFILI = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl__TEM.nc', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl_b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1_TEM.nc']
H0As = ['/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_regrid_0.25x0.25_onpres/*h0a*.nc',\
       '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/atm/hist_regrid_0.25x0.25_onpres/*h0a*.nc']
ALIASES = ['ctrl', 'seed']
DENSFILI = 'tcdens.nc'

PSI_NAME = 'PSI_resid'
pres_name = 'plev'
DIAB_VARS = ['DTCOND', 'QRL', 'QRS'] #TODO: check DTV

LATLAB = linevslat.LATLAB

def main():
   sfs = [xr.open_dataset(os.path.join(DIRIO, fil))[PSI_NAME] for fil in SFFILI]
   pltsf = sfs[1] if DIFF else sfs[0] #streamfuncs are saved as diffs
   expo = 9 if DIFF else 10

   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (8, 12)
   subplot_kw = dict(xlim=(-1, 1))
   fig, axes = plt.subplots(3, 2, sharex=True, subplot_kw=subplot_kw, height_ratios=[2, 1, 1], width_ratios=[15, 1])
   axes[1, 1].axis('off'); axes[2, 1].axis('off')

   contourfkwargs = {'cmap': 'bwr' if DIFF else 'coolwarm', 'levels': 2.**np.arange(-2, 7, 1), 'norm': colors.SymLogNorm(2.**-1)}
   contourfkwargs['levels'] = np.concatenate((-contourfkwargs['levels'][::-1], contourfkwargs['levels']))
   contourkwargs = {'colors': 'black', 'levels': 2.**np.arange(-2, 7, 1)}
   contourkwargs['levels'] = np.concatenate((-contourkwargs['levels'][::-1], contourkwargs['levels']))
   clabelkwargs = {'inline': 1, 'fontsize': 10, 'colors': 'black', 'fmt': '%.1f'}

   ax = axes[0, 0]
   CSF = ax.contourf(np.sin(np.deg2rad(pltsf.lat)), pltsf[pres_name], pltsf.values / 10**expo, **contourfkwargs)
   CS1 = ax.contour(np.sin(np.deg2rad(pltsf.lat)), pltsf[pres_name], pltsf.values / 10**expo, **contourkwargs) 
   ax.set_ylabel('Pressure [hPa]')
   ax.set_yticks(np.arange(100, 1001, 100))
   plt.colorbar(CSF, cax=axes[0, 1])
   ax.set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
   ax.set_yscale('log')
   #ax.set_xlim(np.sin(np.deg2rad(-15)), np.sin(np.deg2rad(90)))
   ax.set_ylim(100, 1000)
   ax.invert_yaxis()
   ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
   ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
   #ax.clabel(CS1, **clabelkwargs)
   ax.set_title('$\\bar{\Psi}^*$ (10$^{%d}$ kg s$^{-1}$)' % expo)
   ax.set_xlabel('Latitude [°]')

   #plt.show()
   plt.savefig(os.path.join(DIRIO, '%s_%s_TEM_ener_tracks.png' % (ALIASES[0], ALIASES[1] if DIFF else '')), bbox_inches='tight')

   buddss = [xr.open_mfdataset(h0).mean(dim=['lon', 'time']) for h0 in H0As] #TODO: weight months by length
   budds = buddss[1] - buddss[0] if DIFF else buddss[0]
   if not DIFF:
      buddss[1].close()

   Tdot3D = sum([budds[dv] for dv in DIAB_VARS])
   contourkwargs = {'colors': 'red', 'levels': np.arange(5e-2, 5.1e-1, 5e-2)}
   contourkwargs['levels'] = np.concatenate((-contourkwargs['levels'][::-1], contourkwargs['levels']))
   CS2 = ax.contour(np.sin(np.deg2rad(Tdot3D.lat)), Tdot3D[pres_name], Tdot3D.values * 86400, **contourkwargs)

   plt.savefig(os.path.join(DIRIO, '%s_%s_TEM_ener_tracks.png' % (ALIASES[0], ALIASES[1] if DIFF else '')), bbox_inches='tight')
   #plt.show()

   integ3d = c.cp / c.g * Tdot3D.integrate(pres_name) #since cumulative not needed trapint(Tdot3D, Tdot3D[pres_name]).isel({pres_name: -1})
   sfcexch = budds['SHFLX'] - budds['LHFLX']
   sumbud = integ3d + sfc

   ax = axes[1, 0]
   ax.plot(np.sin(np.deg2rad(sfcexch.lat)), sfcexch.values, label='sfc')
   ax.plot(np.sin(np.deg2rad(integ3d.lat)), integ3d.values, label='3D cond+lw+sw')
   ax.plot(np.sin(np.deg2rad(sumbud.lat)), sumbud.values, label='sum')
   ax.hlines(0, -1, 1, color='black', linestyle='--')
   ax.legend(framealpha=0.5)
   ax.set_title('Atmo energy budget [W m$^{-2}$]')

   plt.savefig(os.path.join(DIRIO, '%s_%s_TEM_ener_tracks.png' % (ALIASES[0], ALIASES[1] if DIFF else '')), bbox_inches='tight')  
   #plt.show()

   densds = xr.open_dataset(os.path.join(DIRIO, DENSFILI))
   if DIFF:
      densds = densds.sel(run=ALIASES[1]) - densds.sel(run=ALIASES[0])
   else:
      densds = densds.sel(run=ALIASES[0])

   ax = axes[2, 0]
   ax.plot(np.sin(np.deg2rad(densds.lat)), densds.ace.values, label='ACE')
   ax.set_ylabel('ACE')
   ax1 = ax.twinx()
   ax1.plot(np.sin(np.deg2rad(densds.lat)), densds.h6hurr.values, label='6h hurricane fixes', color='orange')
   ax1.set_ylabel('6h hurricane fixes')
   ax1.legend(framealpha=0.5, loc=1)
   ax.set_title(densds.ace.attrs['norm'])
   ax.legend(framealpha=0.5, loc=2)
   ax.hlines(0, -1, 1, color='black', linestyle='--')
   ax1.hlines(0, -1, 1, color='black', linestyle='--')
   #TODO: align zeroes for diff plot

   #plt.show()
   plt.savefig(os.path.join(DIRIO, '%s_%s_TEM_ener_tracks.png' % (ALIASES[0], ALIASES[1] if DIFF else '')), bbox_inches='tight')
   plt.close()

def thta(T, p):
   """
   Returns potential temperature given temp [K] and pressure [hPa].
   """
   return T * (1000 / p)**(c.Rd/c.cp)


if __name__ == '__main__':
   main()
