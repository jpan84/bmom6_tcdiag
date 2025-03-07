import os
import glob
import datetime as dt
import cftime
import numpy as np
import uxarray as ux
import holoviews as hv
from bokeh.models import FixedTicker
import matplotlib.pyplot as plt
import geoviews.feature as gf
import cartopy.crs as ccrs

case = 'b.e23.BMOM.ne120pg3_sx0.66av1.aqua.zmtau.0826'
alias = 'zmtau.0826'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h1i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120pg3_scrip_170417.nc'
#START, END = dt.datetime(5, 2, 15), dt.datetime(5, 3, 1)
START, END = '0005-02-15-21600', '0005-02-28-21600'
LATBNDS = (-30, -10)
LONBNDS = (100, 160)


case = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl2e-5'
alias = 'curl2e-5.1754618'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h1i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
START, END = '0005-02-14-00000', '0005-02-20-00000'
LATBNDS = (-26, -13)
LONBNDS = (-145, -130)
COFVAR = 'PSL'
#COFVAR = 'OMEGA500'
#COFLIM = (-5, 5)
COFLIM = (970, 1020)
#COFLEV = np.arange(-20, 21, 4)
COFLEV = np.concatenante((np.arange(COFLIM[0], 1001, 5), np.arange(1000, 1011, 2), np.arange(1010, 1021)))

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

         da = uxds[COFVAR].sel(time=tt) / 100
         subda = da.subset.bounding_box(LONBNDS, LATBNDS)
         pts = subda.plot.rasterize(title=tstr, height=1200, width=1200)
         pts = pts.opts(cmap='bwr_r', xlim=LONBNDS, ylim=LATBNDS, color_levels=COFLEV, clim=COFLIM, fontsize=hvfont) #symmetric=True
         #uxds.uxgrid.plot()
         hv.save(pts, os.path.join(diro, '%s_%s.png' % (alias, tstr)))


if __name__ == '__main__':
   main()
