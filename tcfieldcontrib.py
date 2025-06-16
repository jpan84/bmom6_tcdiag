#compute the fraction or component of an h1i field that is accounted for by TCs
#Use raw h1i output and NodeFileFilter'ed h1i output

import os
import uxarray as ux
import numpy as np
import sznl_zm_ux
import pltsettings
import matplotlib.pyplot as plt

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1'
DIRO = './tcfieldszm_250416_seed1x1/'
HRAW = 'atm/hist_0010_h1i'
HNFF = 'atm/nff_tcprec'
HTAPE = '*h1i*.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

g = 9.81
WMAT = sznl_zm_ux.wmat
ZMLATS = np.arange(*sznl_zm_ux.zmlats)
SINLAT = np.sin(np.deg2rad(ZMLATS))
LATLAB = sznl_zm_ux.LATLAB
SZNS = sznl_zm_ux.SZNS
lncolors = plt.cm.jet(np.linspace(0, 1, 4))

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)
   pltsettings.set1()

   print('Opening datasets...')
   dsraw = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HRAW, HTAPE))
   dsnff = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HNFF, HTAPE))

   print('\nPRECT...')
   precttcs, precttot = sznl_fields(dsnff, dsraw, 'PRECT')
   print(precttcs.shape)
   plt_sznl(SINLAT, precttot, 'PRECT', '[m s$^{-1}$]', close=False)
   plt_sznl(SINLAT, precttcs, 'PRECT', '[m s$^{-1}$]', linestyle='dashed')
   ###plt_sznl(SINLAT, precttcs / precttot, 'PRECT_TCfrac', '')

   umf = lambda omg: (omg < 0) * -omg / g
   print('\nUMF 500...')
   umf500tcs, umf500tot = sznl_fields(dsnff, dsraw, 'OMEGA500', afunc=umf)
   plt_sznl(SINLAT, umf500tot, 'UMF500', '[kg m$^{-2}$ s$^{-1}$]', close=False)
   plt_sznl(SINLAT, umf500tcs, 'UMF500', '[kg m$^{-2}$ s$^{-1}$]', linestyle='dashed')
   #plt_sznl(SINLAT, umf500tcs / umf500tot, 'UMF500_TCfrac', '')

   print('\nConvective frac 500...')
   _, cf500 = sznl_fields(dsnff, dsraw, 'OMEGA500', afunc=lambda omg: omg < 0)
   plt_sznl(SINLAT, cf500, 'CF500', '')

   print('\nEHF 850...')
   ehf850tcs, ehf850tot = sznl_fields(dsnff, dsraw, 'V850', var2='T850', afunc=lambda x,y: x * y)
   plt_sznl(SINLAT, ehf850tot, 'EHF850', '[K m s$^{-1}$]', close=False)
   plt_sznl(SINLAT, ehf850tcs, 'EHF850', '[K m s$^{-1}$]', linestyle='dashed')

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
   plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
   plt.savefig(os.path.join(DIRO, '%s.%s.png' % (CASE.split('.')[-1], name)), bbox_inches='tight')
   if close:
      plt.close()

if __name__ == '__main__':
   main()
