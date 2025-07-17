#compute the fraction or component of an h1i field that is accounted for by TCs
#Use raw h1i output and NodeFileFilter'ed h1i output

import os
import operator
import uxarray as ux
import numpy as np
from sznl_funcs import monthly2sznl, stack_hemi_sznl
import pltsettings
import matplotlib.pyplot as plt

ARCHV = '/glade/derecho/scratch/jpan/archive/'
ALIAS = '250417_ctrl'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s' % ALIAS
DIRO = './tcfields2mps_%s/' % ALIAS
HRAWS = 'atm/hist_0010_h1i/*h1i*.nc'
HMASK = 'atm/nff_2mps/*.h1i.0010*'
HNORM = 'atm/hist_norms/yhourmean_h1i.nc' #norms e.g., yhourmean
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
MSKNM = 'TC_R2'

print(CASE, DIRO)

zmlats = (-90, 90, 0.5)
LATLAB = np.arange(-90, 91, 30)
g = 9.81

ZMLATS = np.arange(*zmlats)
SINLAT = np.sin(np.deg2rad(ZMLATS))
LSTYS = ['solid', 'dashed']
LCLRS = ['blue', 'orange']

YLIM = {'UMF500_TCfrac': (0, 0.6), 'PRECT_TCfrac': (0, 0.6), 'UMF500': (0, .012), 'PRECT': (0, 2.3e-7), 'CF500': (0, 1)}

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)
   pltsettings.set1()

   print('Opening datasets...')
   dsraws = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HRAWS))
   dsmask = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HMASK))
   dsnorm = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HNORM))

   #TODO: raw vars
   #consider threshold precip rates? (how much of precip at rates heavier than XX do TCs contribute)
   print('Computing precip contrib...')
   rawvars = ['PRECT']
   var1 = 'PRECT'
   flds = [dsraws[var1], dsraws[var1] * dsmask[MSKNM]]
   monzm = [da.groupby('time.month').mean().zonal_mean(lat=zmlats) for da in flds]
   sznzm = [stack_hemi_sznl(monthly2sznl(da)) for da in monzm]
   ratzm = [sznzm[1] / sznzm[0]]

   print('Plotting precip contrib...')
   plt.rcParams['figure.figsize'] = (10, 12)
   subplot_kw = dict(xlim=(-1, 1), sharex=True)
   fig, axes = plt.subplots(2, 1, subplot_kw=subplot_kw)

   plt.suptitle(var1)
   for ii, pltda in enumerate(sznzm):
      for tt, szn in enumerate(pltda['season']):
         axes[0].plot(SINLAT, pltda.sel(season=szn), label=str(szn.values), color=LCLRS[tt], linestyle=LSTYS[ii])
         if not ii:
            axes[1].plot(SINLAT, ratzm, color=LCLRS[tt])
   axes[0].legend()
   axes[0].set_xticks(SINLAT, labels=LATLAB)
   fig.tight_layout()
   plt.savefig(os.path.join(DIRO, '%s.png' % var1))
   plt.close()

   print('tcfieldcontrib.py done.')
   exit()


   #TODO: directional raw vars
   dirvars = ['SHFLX', 'LHFLX', 'TAUX']#, 'UMF500', 'DMF500', 'UMF850', 'DMF850', 'TAUY'] #consider turning zonal stresses into torques
   for dv in dirvars:
      print('Working on directional var %s...' % dv)
      plt.rcParams['figure.figsize'] = (10, 12)
      subplot_kw = dict(xlim=(-1, 1), sharex=True)
      fig, axes = plt.subplots(2, 1, subplot_kw=subplot_kw)
      axes[0].hlines(0, -1, 1, colors='black', linestyles='dotted')
      for jj, op in enumerate([operator.gt, operator.lt]):
         sgned = dsraws[dv] * op(dsraws[dv], 0)
         flds = [sgned, sgned * dsmask[MSKNM]]
         monzm = [da.groupby('time.month').mean().zonal_mean(lat=zmlats) for da in flds]
         sznzm = [stack_hemi_sznl(monthly2sznl(da)) for da in monzm]
         ratzm = [sznzm[1] / sznzm[0]]
         for tt, szn in enumerate(pltda['season']):
            axes[1].plot(SINLAT, ratzm, color=LCLRS[tt], linestyle=LSTYS[jj])
            for ii, pltda in enumerate(sznzm):
               axes[0].plot(SINLAT, pltda.sel(season=szn), label=str(szn.values) if jj else None, color=LCLRS[tt], linestyle=LSTYS[ii])


   #TODO: meridional eddy fluxes

   #TODO: vertical eddy fluxes

   #TODO: Plot total field, TC field, and fraction

   ###print('\nPRECT...')
   ###precttcs, precttot = sznl_fields(dsnff, dsraw, 'PRECT')
   ###print(precttcs.shape)
   ###plt_sznl(SINLAT, precttot, 'PRECT', '[m s$^{-1}$]', close=False)
   ###plt_sznl(SINLAT, precttcs, 'PRECT', '[m s$^{-1}$]', linestyle='dashed')
   ###plt_sznl(SINLAT, precttcs / precttot, 'PRECT_TCfrac', '', linestyle='dotted')

   ###umf = lambda omg: (omg < 0) * -omg / g
   ###print('\nUMF 500...')
   ###umf500tcs, umf500tot = sznl_fields(dsnff, dsraw, 'OMEGA500', afunc=umf)
   ###plt_sznl(SINLAT, umf500tot, 'UMF500', '[kg m$^{-2}$ s$^{-1}$]', close=False)
   ###plt_sznl(SINLAT, umf500tcs, 'UMF500', '[kg m$^{-2}$ s$^{-1}$]', linestyle='dashed')
   ###plt_sznl(SINLAT, umf500tcs / umf500tot, 'UMF500_TCfrac', '', linestyle='dotted')

   print('\nConvective frac 500...')
   _, cf500 = sznl_fields(dsnff, dsraw, 'OMEGA500', afunc=lambda omg: omg < 0)
   plt_sznl(SINLAT, cf500, 'CF500', '')

   #incorrect considered total instead of just eddy flux
   ###print('\nEHF 850...')
   ###ehf850tcs, ehf850tot = sznl_fields(dsnff, dsraw, 'V850', var2='T850', afunc=lambda x,y: x * y)
   ###plt_sznl(SINLAT, ehf850tot, 'EHF850', '[K m s$^{-1}$]', close=False)
   ###plt_sznl(SINLAT, ehf850tcs, 'EHF850', '[K m s$^{-1}$]', linestyle='dashed')

   #TODO: consider eddy moisture flux

#if var2 is (not) None, then afunc must take (2) 1 args.
def sznl_fields(dsnff, dsraw, var1, var2=None, afunc=lambda x: x):
   fldtot, fldtcs = dsraw[var1], dsnff[var1]
   print('Applying func...')
   if not var2 is None:
      fldtot = afunc(fldtot, dsraw[var2])
      fldtcs = afunc(fldtcs, dsnff[var2])
   else:
      fldtot = afunc(fldtot)
      fldtcs = afunc(fldtcs)

   print('Monthly means...')
   totmon = fldtot.groupby('time.month').mean()
   tcsmon = fldtcs.groupby('time.month').mean()
   print('Zonal means...')
   totzm = totmon.zonal_mean(lat=ZMLATS).transpose('month', ...)
   tcszm = tcsmon.zonal_mean(lat=ZMLATS).transpose('month', ...)

   return WMAT @ tcszm.data, WMAT @ totzm.data

#vals should have shape (season, lat)
def plt_sznl(sinlats, arr, name, units, linestyle='solid', close=True):
   for tt in range(arr.shape[0]):
      plt.plot(sinlats, arr[tt, :], label=SZNS[tt], color=lncolors[tt], linestyle=linestyle)
   plt.xlabel('Lat [°]')
   plt.ylabel(name + ' ' + units)
   plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=10))
   plt.xlim(-1, 1)
   plt.ylim(*YLIM[name])
   plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
   plt.savefig(os.path.join(DIRO, '%s.%s.png' % (CASE.split('.')[-1], name)), bbox_inches='tight')
   if close:
      plt.close()

if __name__ == '__main__':
   main()
