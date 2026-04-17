import os
import pandas as pd
import numpy as np
import uxarray as ux
import xarray as xr
import matplotlib.pyplot as plt

RSDIR = '/glade/derecho/scratch/jpan/unseed_restarts_paper1/'
EVNTS = '~/aquaptc/tempest/250415_unseed_production_unseed_events.parquet'
GRIDF = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

DTSTR = '0005-08-07-00000'
CLAT, CLON = 18.00722, 98.79274

MODSTR = '*%s.nc'
ORISTR = '*%s.nc.ORIG.nc'

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
   az_mean = lambda da: da.azimuthal_mean((true_ctr.face_lon, true_ctr.face_lat), 5, 0.25)

   test_p = mirror_azim_mean(az_mean(p_lev).isel(state=0).squeeze())
   print(test_p)
   plt.contour(test_p['radius'], lev_coord.isel(state=0), test_p, levels=np.arange(1.5e4, 9.6e4, 5e3), colors='black')
   plt.ylim(1000, 70)
   plt.yscale('log')
   plt.show()

def mirror_azim_mean(am):
   flip = am.isel(radius=slice(-1, None, -1))
   flip = flip.assign_coords(radius=-flip['radius'])
   return xr.concat([flip, am], dim='radius')

if __name__ == '__main__':
   main()
