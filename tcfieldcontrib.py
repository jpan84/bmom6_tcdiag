#compute the fraction or component of an h1i field that is accounted for by TCs
#Use raw h1i output and NodeFileFilter'ed h1i output

import os
import uxarray as ux
import numpy as np
import sznl_zm_ux

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl'
DIRO = './tcfieldszm_250417_ctrl/'
HRAW = 'atm/hist_0010_h1i'
HNFF = 'atm/nff_tcprec'
HTAPE = '*h1i*.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

g = 9.81
WMAT = sznl_zm_ux.wmat
ZMLATS = sznl_zm_ux.zmlats

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)
   pltsettings.set1()

   dsraw = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HRAW, HTAPE))
   dsnff = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HNFF, HTAPE))

#if var2 is (not) None, then afunc must take (2) 1 args.
def computefields(dsnff, dsraw, var1, var2=None, afunc=lambda x: x):
   fldtot, fldtcs = dsraw[var1], dsnff[var1]
   if not var2 is None:
      fldtot = afunc(fldtot, dsraw[var2])
      fldtcs = afunc(fldtcs, dsnff[var2])
   else:
      fldtot = afunc(fldtot)
      fldtcs = afunc(fldtcs)

   totmon = fldtot.groupby('time.month').mean()
   tcsmon = fldtcs.groupby('time.month').mean()
   totzm = totmon.zonal_mean(lat=ZMLATS).transpose('month', ...)
   tcszm = tcsmon.zonal_mean(lat=ZMLATS).transpose('month', ...)

   return WMAT @ tcszm, WMAT @ totzm

if __name__ == '__main__':
   main()
