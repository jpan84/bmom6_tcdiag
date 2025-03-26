import os
import glob
import datetime as dt
import cftime
import numpy as np
import uxarray as ux
import holoviews as hv
from holoviews.operation import contours as hvcontours
from bokeh.models import FixedTicker

COFVAR = 'PS'
#COFVAR = 'OMEGA500'
#COFLIM = (-5, 5)
#COFLIM = (970, 1030)
#COFLEV = np.concatenate((np.arange(COFLIM[0], 1010, 5), np.arange(1010, 1021, 2), np.arange(1020, 1031)))
COFLIM = (1004, 1024)
COFLEV = np.arange(1004, 1025, 1)
CONVAR = 'VBOT'
CONLEV = np.arange(-20, 0, 4)
CONLEV = np.concatenate((CONLEV, -CONLEV[::-1]))

case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl2e-5'
alias = 'curl2e-5.1754618.h2i'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h2i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
START, END = '0005-02-14-03600', '0005-02-19-03600'
LATBNDS = (-30, -10)
LONBNDS = (-150, -127)

case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl2e-5'
alias = 'curl2e-5.1727914.h2i'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h2i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
START, END = '0005-02-07-03600', '0005-02-10-03600'
LATBNDS = (1, 21)
LONBNDS = (29, 49)

###case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl2e-5'
###alias = 'curl2e-5.1724350'
###diro = 'tcmov_%s' % alias
###diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
###hstr = '*.h1i.%s.nc'
###grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
###START, END = '0005-02-01-21600', '0005-02-07-00000'
###LATBNDS = (-37, -7)
###LONBNDS = (-169, -139)
###
###case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl2e-5'
###alias = 'curl2e-5.1759110'
###diro = 'tcmov_%s' % alias
###diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
###hstr = '*.h1i.%s.nc'
###grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
###START, END = '0005-02-22-43200', '0005-02-27-00000'
###LATBNDS = (-35, -5)
###LONBNDS = (93, 123)
###
###case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl1deg'
###alias = 'curl1deg.1757806'
###diro = 'tcmov_%s' % alias
###diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
###hstr = '*.h1i.%s.nc'
###grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
###START, END = '0005-02-18-00000', '0005-02-24-00000'
###LATBNDS = (-41, -11)
###LONBNDS = (-11, 19)
###
###case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl1deg'
###alias = 'curl1deg.1729575'
###diro = 'tcmov_%s' % alias
###diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
###hstr = '*.h1i.%s.nc'
###grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
###START, END = '0005-02-08-00000', '0005-02-14-00000'
###LATBNDS = (-39, -9)
###LONBNDS = (-93, -63)

case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl1deg'
alias = 'curl1deg.1723746.h2i'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h2i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
START, END = '0005-02-01-03600', '0005-02-05-03600'
LATBNDS = (0, 20)
LONBNDS = (67, 87)

###case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl1deg2e-5'
###alias = 'curl1deg2e-5.NHTC'
###diro = 'tcmov_%s' % alias
###diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
###hstr = '*.h1i.%s.nc'
###grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
###START, END = '0005-02-18-43200', '0005-02-25-00000'
###LATBNDS = (25, 70)
###LONBNDS = (100, 145)
###
###case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl1deg'
###alias = 'curl1deg.1722980'
###diro = 'tcmov_%s' % alias
###diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
###hstr = '*.h1i.%s.nc'
###grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
###START, END = '0005-02-01-21600', '0005-02-04-00000'
###LATBNDS = (-30, -10)
###LONBNDS = (-8, 12)
###
###case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl1deg'
###alias = 'curl1deg.00050202-000000-ncol296899'
###diro = 'tcmov_%s' % alias
###diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
###hstr = '*.h1i.%s.nc'
###grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
###START, END = '0005-02-01-21600', '0005-02-04-00000'
###LATBNDS = (-35, -5)
###LONBNDS = (-180, -155)
###
###case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl1deg'
###alias = 'curl1deg.00050214-000000-ncol430406'
###diro = 'tcmov_%s' % alias
###diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
###hstr = '*.h1i.%s.nc'
###grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
###START, END = '0005-02-11-00000', '0005-02-17-00000'
###LATBNDS = (-28, 2)
###LONBNDS = (-64, -34)

case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250321_unseed_400.15.8e-5'
alias = '250321_unseed.3354390'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h1i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
START, END = '0005-02-12-00000', '0005-02-16-00000'
LATBNDS = (-43, -13)
LONBNDS = (47, 77)

case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250321_unseed_400.15.8e-5'
alias = '250321_unseed.3358807'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h1i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
START, END = '0005-02-22-00000', '0005-02-26-00000'
LATBNDS = (-35, -5)
LONBNDS = (-64, -34)

'''
case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_seed1x1'
alias = 'seed1x1.minp32.73'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h1i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
START, END = '0001-02-01-21600', '0001-02-05-00000'
LATBNDS = (5, 35)
LONBNDS = (-115, -85)

case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_seed1x1'
alias = 'seed1x1.minp31.59'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h1i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
START, END = '0001-02-03-00000', '0001-02-07-00000'
LATBNDS = (-21, 9)
LONBNDS = (152, 180)

case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_seed1x1'
alias = 'seed1x1.minp15.15'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h1i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
START, END = '0001-02-03-00000', '0001-02-07-00000'
LATBNDS = (-10, 20)
LONBNDS = (77, 107)
'''

hvfont = {
    'title': 24, 
    'labels': 24, 
    'xticks': 24, 
    'yticks': 24,
    'ticks': 24}

def main():
   if not os.path.exists(diro):
      os.makedirs(diro)
   mygl = lambda fill: glob.glob(os.path.join(diri, hstr % fill))
   allf = sorted(mygl('*'))
   gls, gle = mygl(START)[0], mygl(END)[0]
   ids, ide = allf.index(gls), allf.index(gle)
   allf = allf[ids: ide + 1]

   for ff in allf:
      uxds = ux.open_dataset(grid, ff)
      for tt in uxds.time:
         dt = cftime.num2date(tt, tt.units, calendar=tt.calendar)
         tstr = dt.strftime('%Y-%m-%d-%H')
         print(tstr)

         framekwargs = dict(height=800, width=800, xlim=LONBNDS, ylim=LATBNDS)

         da = uxds[COFVAR].sel(time=tt) / 100
         dac = uxds[CONVAR].sel(time=tt)
         subda = da.subset.bounding_box(LONBNDS, LATBNDS)
         subdac = dac.subset.bounding_box(LONBNDS, LATBNDS)
         crast = subdac.plot.rasterize(method='polygon', backend='bokeh')
         pts = subda.plot.rasterize(title='%s min=%.1f' % (tstr, subda.min().values), **framekwargs)
         pts = pts.opts(cmap='bwr_r', color_levels=list(COFLEV), clim=COFLIM, fontsize=hvfont) #symmetric=True
         pts *= hvcontours(crast, levels=CONLEV).opts(show_legend=False, cmap='PRGn', color_levels=list(CONLEV), **framekwargs)
         hv.save(pts, os.path.join(diro, '%s_%s.png' % (alias, tstr)))


if __name__ == '__main__':
   main()
