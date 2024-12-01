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

'''
fili = '/glade/derecho/scratch/jpan/archive/QPC4-ne30np4-aquap10-ctrl/atm/hist/QPC4-ne30np4-aquap10-ctrl.cam.h1.1999-10-04-00000.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne30np4_091226_pentagons.nc'
LATBNDS = (15, 25)
LONBNDS = (-145, -135)

fili = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120pg3_sx0.66av1.aqua.reitab.0819/atm/hist/b.e23.BMOM.ne120pg3_sx0.66av1.aqua.reitab.0819.cam.h1i.0004-10-18-21600.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120pg3_scrip_170417.nc'
LATBNDS = (18, 28)
LONBNDS = (-175, -165)
'''
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

case = 'b.e23.BMOM.ne120pg3_sx0.66av1.aqua.reitab.0819'
alias = 'reitab.0819'
diro = 'tcmov_%s' % alias
diri = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/' % case
hstr = '*.h1i.%s.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120pg3_scrip_170417.nc'
START, END = '0004-01-01-21600', '0004-01-31-21600'
LATBNDS = (-50, -20)
LONBNDS = (-60, 60)

COFVAR = 'PSL'
#COFVAR = 'OMEGA500'
COFLIM = (-5, 5)
COFLIM = (950, 1020)
#COFLEV = np.arange(-20, 21, 4)
hvfont = {
    'title': 24, 
    'labels': 24, 
    'xticks': 24, 
    'yticks': 24,
    'ticks': 24}

def main():
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
         #print(uxds.dims)

         #limidx = np.where((uxds.lat >= LATBNDS[0]) & (uxds.lat <= LATBNDS[1]) & (uxds.lon >= LONBNDS[0]) & (uxds.lon <= LONBNDS[1]))[0]

         '''
         #plt.rc('font', size=20)
         #plt.figure(figsize=(12, 8))
         rast = uxds[COFVAR].sel(time=tt).plot.polygons(backend='bokeh')
         #rast = rast.opts(cmap='bwr', symmetric=True, xlim=LONBNDS, ylim=LATBNDS, color_levels=10, clim=COFLIM,
         #   line_width=0.1)# * features
         rast = rast.opts(cmap='bwr', symmetric=True, xlim=LONBNDS, ylim=LATBNDS, color_levels=10, clim=COFLIM,
             frame_width=1200, frame_height=800)# * features
         hv.save(rast, os.path.join(diro, '%s_%s.png' % (alias, tstr)))
         #plt.close()
         '''

         pts = (uxds[COFVAR].sel(time=tt) / 100).plot.points(title=tstr, height=600, width=1600, size=2)
         pts = pts.opts(cmap='bwr_r', xlim=LONBNDS, ylim=LATBNDS, color_levels=7, clim=COFLIM, fontsize=hvfont) #symmetric=True
         #uxds.uxgrid.plot()
         hv.save(pts, os.path.join(diro, '%s_%s.png' % (alias, tstr)))


if __name__ == '__main__':
   main()
