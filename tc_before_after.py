import os
import re
import numpy as np
import uxarray as ux
import holoviews as hv
from bokeh.models import FixedTicker
import geoviews.feature as gf
import cartopy.crs as ccrs


BBOX_DEG = 10.
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all'
OUTDIR = 'bef_aft.250129_unseed_all/'
seedlog = open('/glade/u/home/jpan/work/MOM6_CASEDIRS/250129_ne120np4_unseed_all.out4', 'r')
RESTDIR = '/glade/derecho/scratch/jpan/archive/%s/rest/%s/'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
ps_pattern = r"^\[(1\d{5}|9\d{4})\."

def main():
   if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

   mkplt = False
   clat, clon, psmin, tid = None, None, None, None
   dtstr, orinc, modnc = None, None, None
   for ln in seedlog:
      spl = ln.strip('\n').split(' ')
      if not len(spl[0]):
         continue
      if spl[0] == 'sedding':
         #211 sedding lat/lon: -8.347259 48.750000
         clon = float(spl[-1])
         clat = float(spl[-2])
      elif spl[:2] == ['Debugging', 'restart']:
         #212 Debugging restart file /glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc at time 5407747
         tid = int(spl[-1])
      elif spl[0][0] == "'":
         #200 '/glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc' -> '/glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc.ORIG.nc'
         dtstr = re.search(r"(\d{4}-\d{2}-\d{2}-\d{5})", spl[0]).group()
         print(dtstr)
         modnc = os.path.basename(spl[0].strip("'"))
         orinc = os.path.basename(spl[-1].strip("'"))
         #print(modnc, orinc)
      elif spl[0][0] == '[':
         if re.match(ps_pattern, spl[0]):
            psmin = float(spl[0].strip('[')) / 100
            mkplt = True #make plot now that we have all the info for 1 storm

      if mkplt:
         print('Plotting %d...' % tid)
         rdir = RESTDIR % (CASE, dtstr)
         ori = os.path.join(rdir, orinc)
         mod = os.path.join(rdir, modnc)
         orids = ux.open_dataset(grid, ori)
         modds = ux.open_dataset(grid, mod)
         print('\topened')

         latbnds = [clat - BBOX_DEG, clat + BBOX_DEG]
         lonbnds = [clon - BBOX_DEG, clon + BBOX_DEG]
         lonbnds = [(lon + 180) % 360 - 180 for lon in lonbnds]
         extrps = lambda ds: ds['PSDRY'].isel(time=0) / 100
         selbbox = lambda da: da.subset.bounding_box(lonbnds, latbnds)
         print('\tselected')

         orivar = selbbox(extrps(orids))
         modvar = selbbox(extrps(modds))
         difvar = modvar - orivar
         print('\tdiffed')

         rast = difvar.plot.rasterize(method='polygon', backend='bokeh')#, projection=ccrs.PlateCarree())
         #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
         #print(type(rast))
         rast = rast.opts(cmap='bwr', xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400, symmetric=True)
         hv.save(rast, os.path.join(OUTDIR, 'difps_ux_%s_%s.png' % (tid, dtstr)))

         clim = (psmin // 5 * 5, 1020)
         rast = orivar.plot.rasterize(method='polygon', backend='bokeh')#, projection=ccrs.PlateCarree())
         #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
         rast = rast.opts(cmap='BuPu', clim=clim, xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400)
         hv.save(rast, os.path.join(OUTDIR, 'orips_ux_%s_%s.png' % (tid, dtstr)))

         rast = modvar.plot.rasterize(method='polygon', backend='bokeh')#, projection=ccrs.PlateCarree())
         #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
         rast = rast.opts(cmap='BuPu', clim=clim, xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400)
         hv.save(rast, os.path.join(OUTDIR, 'modps_ux_%s_%s.png' % (tid, dtstr)))

         mkplt = False

   exit()


   LATBNDS = (-18, 2)
   LONBNDS = (-8.5, 11.5)
   
   origds = ux.open_dataset(grid, os.path.join(RESTDIR, orig))
   modds = ux.open_dataset(grid, os.path.join(RESTDIR, mod))
   print(list(modds.data_vars))
   #print(uxds.OMEGA500)
   
   
   #https://uxarray.readthedocs.io/en/latest/user-guide/subset.html
   bbox_subset_nodes = (origds['PSDRY'].isel(time=0) / 100).subset.bounding_box(
       LONBNDS,
       LATBNDS,
   )
   bbox_subset_nodes.attrs['units'] = 'hPa'
   
   rast = bbox_subset_nodes.plot.rasterize(method='polygon', backend='bokeh')#, projection=ccrs.PlateCarree())
   #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
   print(type(rast))
   rast = rast.opts(cmap='BuPu', color_levels=list(np.arange(930, 1030, 10)), xlim=LONBNDS, ylim=LATBNDS, aspect='square', frame_width=400)#symmetric=True 
   hv.save(rast, 'subset_tcsnap_ux.png')

if __name__ == '__main__':
   main()
