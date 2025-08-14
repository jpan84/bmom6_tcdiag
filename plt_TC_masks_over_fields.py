import uxarray as ux
import numpy as np
import pandas as pd
import os
import glob
import cftime

import holoviews as hv
from holoviews.operation import contours as hvcontours

CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl'
ATM = '/glade/derecho/scratch/jpan/archive/%s/atm'
H1I = os.path.join(ATM % CASE, 'hist/') #h1i path
MASKS = os.path.join(ATM % CASE, 'nff_2mps/') #masks path
PARQ = '/glade/u/home/jpan/aquaptc/tempest/250417_ctrl.parquet' #TC parquet
CAMGRID = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
DIRO = './TC_mask_plots_R2/'

STRF = lambda dtobj: f'*{dtobj.year:04}-{dtobj.month:02}-{dtobj.day:02}-{3600*dtobj.hour:05}*'

BBDEG = 12.
MNM = 'TC_R2'

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   tcdf = pd.read_parquet(PARQ)
   print(tcdf)

   for ii, row in tcdf.iterrows(): #For timestep in TC parquet

      if np.isnan(row['year']) or row['year'] != 2010.:
         continue
      truedt = cftime.DatetimeNoLeap(int(row['year']) - 2000, int(row['month']), int(row['day']), hour=int(row['hour']))
      print('\n', truedt, STRF(truedt))

      #open h1i
      hglob = glob.glob(os.path.join(H1I, STRF(truedt)))
      print(hglob)
      #open mask
      mglob = glob.glob(os.path.join(MASKS, STRF(truedt)))
      print(mglob) 
      if len(hglob) == 0 or len(mglob) == 0:
         continue
      hds = ux.open_dataset(CAMGRID, hglob[0])
      mds = ux.open_dataset(CAMGRID, mglob[0])

      #print(hds)
      #print(mds)
      print(row['stmnum'], row['lon'], row['lat'])
      lonshift = lambda lon: (lon + 180) % 360 - 180
      clon = row['lon']
      lonbnds = (lonshift(clon - BBDEG), lonshift(clon + BBDEG))
      latbnds = (row['lat'] - BBDEG, row['lat'] + BBDEG)
      msub = mds[MNM].subset.bounding_box(lonbnds, latbnds).squeeze()
      print(msub)

      framekwargs = dict(height=800, width=800, xlim=lonbnds, ylim=latbnds)
      #mrast = msub.plot.rasterize(method='polygon', backend='bokeh')
      mrast = msub.plot.polygons(rasterize=True, backend='bokeh')
      mcont = hvcontours(mrast).opts(show_legend=False, **framekwargs)
      hv.save(mcont, os.path.join(DIRO, '%s_stm%d_%s.png' % (STRF(truedt)[1:-1], row['stmnum'], MNM)))
      #For field in ['U850', 'V850', 'PRECT', 'OMEGA500', 'PS']
      for fld in ['U850', 'V850', 'PRECT', 'OMEGA500', 'PS']:
         print('\t', fld)
         #contourf/rasterize plot field
         hsub = hds[fld].subset.bounding_box(lonbnds, latbnds).squeeze()
         hrast = hsub.plot.polygons(rasterize=True, backend='bokeh', **framekwargs)
         if fld == 'OMEGA500':
            hrast = hrast.opts(clim=(-hsub.max().values.item(), hsub.max().values.item()))
         if fld == 'PRECT':
            hrast = hrast.opts(cnorm='log', clim=(1e-8, 1e-5))
         if fld in ['PRECT', 'PS']:
            hrast = hrast.opts(cmap='viridis')
         if fld in ['OMEGA500', 'U850', 'V850']:
            hrast = hrast.opts(symmetric=True, cmap='seismic')
         #contour mask
         hrast *= mcont
         hv.save(hrast, os.path.join(DIRO, '%s_stm%d_%s.png' % (STRF(truedt)[1:-1], row['stmnum'], fld)))

if __name__ == '__main__':
   main()
