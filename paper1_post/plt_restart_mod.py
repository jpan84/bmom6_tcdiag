import os
import pandas as pd
import numpy as np
import uxarray as ux
import xarray as xr
import matplotlib.pyplot as plt

RSDIR = '/glade/derecho/scratch/jpan/unseed_restarts_paper1/'
EVNTS = '~/aquaptc/tempest/250415_unseed_production_unseed_events.parquet'
GRIDF = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
GRIDDELTA = 0.25

DTSTR = '0005-08-07-00000'
CLAT, CLON = 18.00722, 98.79274
RADO = 2. #gcd

#DTSTR = '0005-09-03-00000'
#CLAT, CLON = 29.45009, 141.5423
#RADO = 3. #gcd

MODSTR = '*%s.nc'
ORISTR = '*%s.nc.ORIG.nc'

PLEVS = levels=np.arange(-1e4, 0, 5e2)
PLEVS = np.concatenate((PLEVS, -PLEVS[::-1]))

TLEVS = levels=np.arange(-5, 0, .5)
TLEVS = np.concatenate((TLEVS, -TLEVS[::-1]))

def main():
   orids = ux.open_mfdataset(GRIDF, os.path.join(RSDIR, ORISTR % DTSTR)).expand_dims(state=['before'])
   modds = ux.open_mfdataset(GRIDF, os.path.join(RSDIR, MODSTR % DTSTR)).expand_dims(state=['after'])
   modds.uxgrid = orids.uxgrid
   ds = ux.concat([orids, modds], dim='state')

   ps = ds['PSDRY'] + ds['dpQ'].sum('lev')
   aiterm = ds['hyai'] * ds['P0']
   biterm = ds['hybi'] * ps
   p_ilev = aiterm + biterm
   dp3d = p_ilev.diff('ilev').rename(dict(ilev='lev')).assign_coords(lev=ds['lev'])
   q = ds['dpQ'] / dp3d
   amterm = ds['hyam'] * ds['P0']
   bmterm = ds['hybm'] * ps
   p_lev = amterm + bmterm
   lev_coord = 1000. * (ds['hyam'] + ds['hybm'])
   ds = ds.assign(q=q, p_lev=p_lev)
   #print(p_ilev / ps)

   #TODO: compute tangential wind

   true_ctr = ds.uxgrid.subset.nearest_neighbor((CLON, CLAT), 1)
   az_mean = lambda da: da.azimuthal_mean((true_ctr.face_lon, true_ctr.face_lat), RADO, GRIDDELTA)

   #test_p = mirror_azim_mean(az_mean(p_lev).isel(state=0).squeeze())
   #test_p -= test_p.isel(radius=-1)
   #test_vt = mirror_azim_mean(az_mean(uv_to_tang(ds, 'U', 'V', (CLON, CLAT), radbound=RADO)).isel(state=0).squeeze())
   #print(test_vt)
   ##plt.contour(test_p['radius'], lev_coord.isel(state=0), test_p, levels=np.arange(1.5e4, 9.6e4, 5e3), colors='black')
   #plt.contourf(test_vt['radius'], lev_coord.isel(state=0), test_vt, levels=np.arange(-50, 51, 5), cmap='PRGn')
   #plt.colorbar()
   #plt.contour(test_p['radius'], lev_coord.isel(state=0), test_p, levels=PLEVS, colors='black')
   #plt.ylim(1000, 70)
   #plt.yscale('log')
   #plt.show()

   xsect = lambda da: da.cross_section(start=(CLON - RADO, CLAT), end=(CLON + RADO, CLAT), steps=int(2*RADO/GRIDDELTA))

   psec = xsect(p_lev.isel(state=0).squeeze()) - az_mean(p_lev).isel(state=0).squeeze().isel(radius=-1).data[:, None]
   vsec = xsect(ds['V'].isel(state=0).squeeze())
   plt.contourf(vsec['lon'], lev_coord.isel(state=0), vsec, levels=np.arange(-50, 51, 5), cmap='PRGn')
   plt.colorbar()
   plt.contour(psec['lon'], lev_coord.isel(state=0), psec, levels=PLEVS, colors='black')
   plt.ylim(1000, 70)
   plt.yscale('log')
   plt.show()

   tsec = xsect(ds['T'].isel(state=0).squeeze()) - az_mean(ds['T']).isel(state=0).squeeze().isel(radius=-1).data[:, None]
   qsec = xsect(q.isel(state=0).squeeze())
   plt.contourf(qsec['lon'], lev_coord.isel(state=0), qsec, levels=np.arange(2e-3, 3.1e-2, 2e-3), cmap='YlGnBu')
   plt.colorbar()
   plt.contour(tsec['lon'], lev_coord.isel(state=0), tsec, levels=TLEVS, colors='black')
   plt.ylim(1000, 70)
   plt.yscale('log')
   plt.show()


def mirror_azim_mean(am):
   flip = am.isel(radius=slice(-1, None, -1))
   flip = flip.assign_coords(radius=-flip['radius'])
   return xr.concat([flip, am], dim='radius')

def uv_to_tang(uxds, unm, vnm, ctrcoord, radbound=None):
   myu, myv, mylat, mylon = uxds[unm], uxds[vnm], uxds['lat'], uxds['lon']
   if radbound is not None:
      bc = lambda da: da.subset.bounding_circle(ctrcoord, radbound + 1, element='face centers')
      myu, myv, mylat, mylon = bc(myu), bc(myv), bc(mylat), bc(mylon)
   clon, clat = ctrcoord
   d1 = np.sin(np.deg2rad(clat)) * np.cos(np.deg2rad(mylat)) - np.cos(np.deg2rad(clat)) * np.sin(np.deg2rad(mylat)) * np.cos(np.deg2rad(mylon) - np.deg2rad(clon))
   d2 = np.cos(np.deg2rad(clat)) * np.sin(np.deg2rad(mylon) - np.deg2rad(clon))
   d = np.maximum(1e-25, np.sqrt(d1 ** 2. + d2 ** 2.))

   ufac = d1 / d
   vfac = d2 / d

   #print(type(ufac))
   #return ux.UxDataArray(ufac**2+vfac**2, uxgrid=mylat.uxgrid)
   #return -myu * ufac
   return myv * vfac + myu * ufac


if __name__ == '__main__':
   main()
