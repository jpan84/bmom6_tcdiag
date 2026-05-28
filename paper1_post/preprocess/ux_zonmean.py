import re
import os
import sys
#sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
#from paths import CAMGR
import uxarray as ux
import xarray as xr
from dask.diagnostics import ProgressBar

print(sys.argv)
PTHPTR = sys.argv[1]#.replace('\\', '')
HPTH = ''
with open(PTHPTR, 'r') as f:
   HPTH = f.read().strip()
DIRI = os.path.dirname(HPTH)
TAPE = re.findall(r'h[0-9][ia]', HPTH)
#TAPE = '*h0a.00[0-9]*.nc' #file names to process

UGRD = sys.argv[2]
CONS = (sys.argv[3] == 'True') #conservative or not
VARS = sys.argv[4] #example of comma-separated (no space): PS,PRECT,QFLX

LATS = (-90, 90, 0.25)

ds = ux.open_mfdataset(UGRD, HPTH) #os.path.join(DIRI, TAPE))
outds = None
myvars = [str(dv) for dv in ds.data_vars]
if VARS != 'all':
   myvars = VARS.split(',')
for dv in myvars:
   zm = ds[dv].zonal_mean(lat=LATS, conservative=CONS)
   del zm.attrs['zonal_mean'], zm.attrs['conservative']
   if outds is None:
      outds = xr.Dataset(data_vars={dv: zm})
   else:
      outds = outds.assign(variables={dv: zm})

outds = outds.assign_attrs(script_from=sys.argv[0], conservative='True' if CONS else 'False', zmlats=str(LATS), ugrd=UGRD, infiles=TAPE)#, zonal_mean=str(outds.attrs['zonal_mean']))
nameflags = ['uxzm', TAPE, ('' if CONS else 'non') + 'cons'] + [str(itm) for itm in LATS] + ([] if VARS == 'all' else myvars)
with ProgressBar():
   outds.to_netcdf(os.path.join(DIRI, '..', '_'.join(nameflags) + '.nc'))
