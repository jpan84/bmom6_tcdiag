import re
import os
import sys
import uxarray as ux
import xarray as xr
from dask.diagnostics import ProgressBar

print(sys.argv)
PTHPTR = sys.argv[1] #file containing the glob path of model output files
HPTH = ''
with open(PTHPTR, 'r') as f:
   HPTH = f.read().strip()
print(HPTH)
os.remove(PTHPTR)
DIRI = os.path.dirname(HPTH)
TAPE = re.findall(r'h[0-9][ia]', HPTH)[0]

UGRD = sys.argv[2]
CONS = (sys.argv[3] == 'True') #conservative or not
VARS = sys.argv[4] #example of comma-separated (no space): PS,PRECT,QFLX
LATS = tuple([float(el) for el in sys.argv[5].split(':')])

ds = ux.open_mfdataset(UGRD, HPTH)
outds = None
myvars = [str(dv) for dv in ds.data_vars]
if VARS != 'all':
   myvars = VARS.split(',')
for dv in myvars:
   zm = ds[dv].zonal_mean(lat=LATS, conservative=CONS)
   del zm.attrs['zonal_mean'], zm.attrs['conservative']
   oricoords = dict(ds[dv].coords)
   zm = zm.assign_coords(coords=oricoords)

   if outds is None:
      outds = xr.Dataset()
   outds[dv] = zm

outds = outds.assign_attrs(script_from=sys.argv[0], conservative='True' if CONS else 'False', zmlats=str(LATS), ugrd=UGRD, infiles=HPTH)#, zonal_mean=str(outds.attrs['zonal_mean']))
nameflags = ['uxzm', TAPE, ('' if CONS else 'non') + 'cons'] + [str(itm) for itm in LATS] + ([] if VARS == 'all' else myvars)
with ProgressBar():
   outds.to_netcdf(os.path.join(DIRI, '..', '_'.join(nameflags) + '.nc'))
