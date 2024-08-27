import uxarray as ux
import holoviews as hv
from bokeh.models import FixedTicker
import geoviews.feature as gf
import cartopy.crs as ccrs

fili = '/glade/derecho/scratch/jpan/archive/QPC4-ne30np4-aquap10-ctrl/atm/hist/QPC4-ne30np4-aquap10-ctrl.cam.h1.1999-10-04-00000.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne30np4_091226_pentagons.nc'
LATBNDS = (15, 25)
LONBNDS = (-145, -135)

fili = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120pg3_sx0.66av1.aqua.reitab.0819/atm/hist/b.e23.BMOM.ne120pg3_sx0.66av1.aqua.reitab.0819.cam.h1i.0004-10-18-21600.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120pg3_scrip_170417.nc'
LATBNDS = (18, 28)
LONBNDS = (-175, -165)

uxds = ux.open_dataset(grid, fili)
print(list(uxds.data_vars))
#print(uxds.OMEGA500)

rast = uxds.OMEGA500.isel(time=0).plot.rasterize(method='polygon', backend='bokeh')#, projection=ccrs.PlateCarree())
#features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
print(type(rast))
rast = rast.opts(cmap='bwr', symmetric=True, xlim=LONBNDS, ylim=LATBNDS, aspect='square', frame_width=400)# * features
hv.save(rast, 'test_tcsnap_ux_4.png')
