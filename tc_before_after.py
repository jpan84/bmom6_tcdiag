import os
import uxarray as ux
import holoviews as hv
from bokeh.models import FixedTicker
import geoviews.feature as gf
import cartopy.crs as ccrs

RESTDIR = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/rest/0001-02-07-00000/'
orig = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-07-00000.nc.ORIG.nc'
mod = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-07-00000.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

LATBNDS = (-18, 2)
LONBNDS = (-8.5, 11.5)

origds = ux.open_dataset(grid, os.path.join(RESTDIR, orig))
modds = ux.open_dataset(grid, os.path.join(RESTDIR, mod))
print(list(modds.data_vars))
#print(uxds.OMEGA500)


#https://uxarray.readthedocs.io/en/latest/user-guide/subset.html
bbox_subset_nodes = origds['PSDRY'].isel(time=0).subset.bounding_box(
    LONBNDS,
    LATBNDS,
)


rast = bbox_subset_nodes.plot.rasterize(method='polygon', backend='bokeh')#, projection=ccrs.PlateCarree())
#features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
print(type(rast))
rast = rast.opts(cmap='bwr', xlim=LONBNDS, ylim=LATBNDS, aspect='square', frame_width=400)#symmetric=True 
hv.save(rast, 'subset_tcsnap_ux.png')
