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
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250219_unseed_all'
OUTDIR = 'bef_aft.250219_unseed_all/'
seedlog = open('/glade/u/home/jpan/work/MOM6_CASEDIRS/250219_ne120np4_unseed_all.out', 'r')
RESTDIR = '/glade/derecho/scratch/jpan/archive/%s/rest/%s/'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
ps_pattern = r"^\[(1\d{5}|9\d{4})\."

dpmin = 0 #3 #hPa
dTmax = 0 #-1 #K
CEN_DEG = 3 #find largest dp and dT within this many GCD of clon,clat

Uclim = 20

lev = None

def main():
   print('\n\n\nStarting with dp, dT thresholds\n', dpmin, dTmax)
   global lev
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
         #difds = modds - orids
         print('\topened')

         if lev is None:
            lev = modds.hyam * modds.P0 + modds.hybm * modds['PSDRY'].mean()
            lev /= 100.

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

         levi = -1
         oriu = selbbox(orids['U'].isel(lev=levi).isel(time=0))
         oriu = oriu - oriu.mean()
         oriv = selbbox(orids['V'].isel(lev=levi).isel(time=0))
         oriv = oriv - oriv.mean()
         modu = selbbox(modds['U'].isel(lev=levi).isel(time=0))
         modu = modu - oriu.mean()
         modv = selbbox(modds['V'].isel(lev=levi).isel(time=0))
         modv = modv - oriv.mean()

         
         levT = -9
         oriT = selbbox(modds['T'].isel(lev=levT).isel(time=0))

         #TODO: maybe generalize to allow plotting other vars
         dp = float(selcir(difvar).max().values)
         dTarr = selcir(modds.T) - selcir(orids.T)
         dT = float(dTarr.min().values)

         if dp < dpmin or dT > dTmax:
            print('dT or dp threshold not met. Skipping...\n')
            mkplt = False
            continue

         #TODO: fix holoviews unable to plot xlim1 > xlim2 across the antimeridian
         rastkwargs = dict(method='polygon', backend='bokeh')
         framekwargs = dict(xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400)
         plevs = np.arange(psmin // 2 * 2, 1042, 2)
         panels = [[], []]

         panels[0].append(difvar.plot.rasterize(**rastkwargs))
         panels[0][-1] = panels[0][-1].opts(title='dp=%.2f' % dp, cmap='bwr', clim=(-dp, dp), symmetric=True, **framekwargs)
         orip = orivar.plot.rasterize(**rastkwargs)
         modp = modvar.plot.rasterize(**rastkwargs)

         panels[0].append(oriu.plot.rasterize(**rastkwargs))
         panels[0][-1] = panels[0][-1].opts(title='before u, lev %d' % levi, cmap='BrBG', clim=(-Uclim, Uclim), symmetric=True, **framekwargs)
         panels[0][-1] *= hvcontours(orip).opts(show_legend=False, **framekwargs)

         panels[0].append(oriv.plot.rasterize(**rastkwargs))
         panels[0][-1] = panels[0][-1].opts(title='before v', cmap='BrBG', clim=(-Uclim, Uclim), symmetric=True, **framekwargs)
         panels[0][-1] *= hvcontours(orip).opts(show_legend=False, **framekwargs)

         panels[1].append(oriT.plot.rasterize(**rastkwargs))
         panels[1][-1] = panels[1][-1].opts(title='before T, lev=%.2f' % lev[levT], cmap='jet', **framekwargs)
         
         panels[1].append(modu.plot.rasterize(**rastkwargs))
         panels[1][-1] = panels[1][-1].opts(title='after u', cmap='BrBG', clim=(-Uclim, Uclim), symmetric=True, **framekwargs)
         panels[1][-1] *= hvcontours(modp).opts(show_legend=False, **framekwargs)

         panels[1].append(modv.plot.rasterize(**rastkwargs))
         panels[1][-1] = panels[1][-1].opts(title='after v', cmap='BrBG', clim=(-Uclim, Uclim), symmetric=True, **framekwargs)
         panels[1][-1] *= hvcontours(modp).opts(show_legend=False, **framekwargs)

         layout = hv.Layout(panels[0] + panels[1]).cols(3)

         '''
         prasts, urasts, vrasts = [], [], []
         prasts.append(difvar.plot.rasterize(method='polygon', backend='bokeh'))#, projection=ccrs.PlateCarree())
         #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
         #print(type(rast))
         prasts[-1] = prasts[-1].opts(cmap='bwr', clim=(-dp, dp), xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400, symmetric=True)
         #hv.save(prasts[-1], os.path.join(OUTDIR, 'difps_ux_%s_%s.png' % (tid, dtstr)))

         clim = (psmin // 5 * 5, 1020)
         plevs = np.arange(psmin // 2 * 2, 1042, 2)
         prasts.append(orivar.plot.rasterize(method='polygon', backend='bokeh'))#, projection=ccrs.PlateCarree())
         #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
         prasts[-1] = prasts[-1].opts(cmap='BuPu', clim=clim, xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400)
         urasts.append(oriu.plot.rasterize(method='polygon', backend='bokeh'))
         urasts[-1] = urasts[-1].opts(cmap='BrBG', clim=(-20, 20), xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400, symmetric=True)
         vrasts.append(oriv.plot.rasterize(method='polygon', backend='bokeh'))
         vrasts[-1] = vrasts[-1].opts(cmap='PuOr', clim=(-20, 20), xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400, symmetric=True, alpha=0.6)
         #hv.save(prasts[-1], os.path.join(OUTDIR, 'orips_ux_%s_%s.png' % (tid, dtstr)))

         prasts.append(modvar.plot.rasterize(method='polygon', backend='bokeh'))#, projection=ccrs.PlateCarree())
         #features = gf.coastline(projection=ccrs.PlateCarree(), line_width=1, scale='50m')
         prasts[-1] = prasts[-1].opts(cmap='BuPu', clim=clim, xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400)
         urasts.append(modu.plot.rasterize(method='polygon', backend='bokeh'))
         urasts[-1] = urasts[-1].opts(cmap='BrBG', clim=(-20, 20), xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400, symmetric=True)
         vrasts.append(modv.plot.rasterize(method='polygon', backend='bokeh'))
         vrasts[-1] = vrasts[-1].opts(cmap='PuOr', clim=(-20, 20), xlim=tuple(lonbnds), ylim=tuple(latbnds), aspect='square', frame_width=400, symmetric=True, alpha=0.6)
         #hv.save(prasts[-1], os.path.join(OUTDIR, 'modps_ux_%s_%s.png' % (tid, dtstr)))

         #!TEST-CONTOURS
         cntrs = [hvcontours(rast).opts(color='black', show_legend=False, aspect='square', frame_width=400) for rast in prasts]
         cntrs[1] = urasts[0] * vrasts[0] * cntrs[1]
         cntrs[2] = urasts[1] * vrasts[1] * cntrs[2]

         layout = hv.Layout([prasts[0], cntrs[1], cntrs[2]])#hv.Layout([el[0] * el[1] for el in zip(prasts, cntrs)])
         '''
         hv.save(layout, os.path.join(OUTDIR, 'layout_ps_ux_%s_%s.png' % (tid, dtstr)))

         mkplt = False


if __name__ == '__main__':
   main()
