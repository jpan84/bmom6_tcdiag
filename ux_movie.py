import os
import glob
import datetime as dt
import cftime
import uxarray as ux
import holoviews as hv
from bokeh.models import FixedTicker
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

CONVAR = 'PSL'
COFVAR = 'OMEGA500'
COFLIM = (-20, 20)

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
         #print(tstr)

         rast = uxds[COFVAR].isel(time=0).plot.polygons(backend='bokeh')
         rast = rast.opts(cmap='bwr', symmetric=True, xlim=LONBNDS, ylim=LATBNDS, frame_width=1200, frame_height=800)# * features
         hv.save(rast, os.path.join(diro, '%s_%s.png' % (alias, tstr)))

   exit()
   
   uxds = ux.open_dataset(grid, fili)
   print(list(uxds.data_vars))
   #print(uxds.OMEGA500)
   
   rast = uxds.OMEGA500.isel(time=0).plot.polygons(backend='bokeh')#, projection=ccrs.PlateCarree())
   #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
   print(type(rast))
   rast = rast.opts(cmap='bwr', symmetric=True, xlim=LONBNDS, ylim=LATBNDS, aspect='square', frame_width=400)# * features
   hv.save(rast, 'test_tcsnap_ux_5.png')

if __name__ == '__main__':
   main()
