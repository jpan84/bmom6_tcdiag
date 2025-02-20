import os
import re
import numpy as np
import uxarray as ux

import holoviews as hv
from holoviews.operation import contours as hvcontours
from bokeh.models import FixedTicker
import geoviews.feature as gf
import cartopy.crs as ccrs


BBOX_DEG = 10.
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250218_unseed_all'
OUTDIR = 'bef_aft.250218_unseed_all.contours/'
seedlog = open('/glade/u/home/jpan/work/MOM6_CASEDIRS/250218_ne120np4_unseed_all.out2', 'r')
RESTDIR = '/glade/derecho/scratch/jpan/archive/%s/rest/%s/'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
ps_pattern = r"^\[(1\d{5}|9\d{4})\."

dpmin = 0 #3 #hPa
dTmax = 0 #-1 #K
CEN_DEG = 3 #find largest dp and dT within this many GCD of clon,clat

def main():
   print('\n\n\nStarting with dp, dT thresholds\n', dpmin, dTmax)
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
         selcir = lambda da: da.subset.bounding_circle((clon, clat), CEN_DEG)

         orivar = selbbox(extrps(orids))
         modvar = selbbox(extrps(modds))
         print('\tselected')
         difvar = modvar - orivar
         print('\tdiffed')

         #TODO: maybe generalize to allow plotting other vars
         dp = float(selcir(difvar).max().values)
         dTarr = selcir(modds.T) - selcir(orids.T)
         dT = float(dTarr.min().values)

         if dp < dpmin or dT > dTmax:
            print('dT or dp threshold not met. Skipping...\n')
            mkplt = False
            continue

         #TODO: fix holoviews unable to plot xlim1 > xlim2 across the antimeridian
         rasts = []
         rasts.append(difvar.plot.rasterize(method='polygon', backend='bokeh'))#, projection=ccrs.PlateCarree())
         #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
         #print(type(rast))
         rasts[-1] = rasts[-1].opts(cmap='bwr', clim=(-dp, dp), xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400, symmetric=True)
         hv.save(rasts[-1], os.path.join(OUTDIR, 'difps_ux_%s_%s.png' % (tid, dtstr)))

         clim = (psmin // 5 * 5, 1020)
         rasts.append(orivar.plot.rasterize(method='polygon', backend='bokeh'))#, projection=ccrs.PlateCarree())
         #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
         rasts[-1] = rasts[-1].opts(cmap='BuPu', clim=clim, xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400)
         hv.save(rasts[-1], os.path.join(OUTDIR, 'orips_ux_%s_%s.png' % (tid, dtstr)))

         rasts.append(modvar.plot.rasterize(method='polygon', backend='bokeh'))#, projection=ccrs.PlateCarree())
         #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
         rasts[-1] = rasts[-1].opts(cmap='BuPu', clim=clim, xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400)
         hv.save(rasts[-1], os.path.join(OUTDIR, 'modps_ux_%s_%s.png' % (tid, dtstr)))

         #!TEST-CONTOURS
         cntrs = [hvcontours(rast).opts(color='black', show_legend=False) for rast in rasts]

         layout = hv.Layout(cntrs)#hv.Layout([el[0] * el[1] for el in zip(rasts, cntrs)])
         hv.save(layout, os.path.join(OUTDIR, 'layout_ps_ux_%s_%s.png' % (tid, dtstr)))

         mkplt = False


if __name__ == '__main__':
   main()
