import argparse
import os
import sys

# Betacast modules
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'py_functions'))
if module_path not in sys.path:
    sys.path.append(module_path)
from py_seedfuncs import replace_or_add_variable

# Set up argument parsing
parser = argparse.ArgumentParser(description="Read SST from CAM and output max lat to vortex namelist.")
parser.add_argument('--camhi', required=True, help='Path to cam instantaneous history')
parser.add_argument('--vortex_namelist', required=True, help='Path to the namelist file')
args = parser.parse_args()

ptho = args.vortex_namelist
camhi = args.camhi

import numpy as np
import uxarray as ux

sstlats = (-45, 45, 0.5)
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne30np4_091226_pentagons.nc'

def main():
   print('Finding SST max lat in', sstlats)
   ts = ux.open_dataset(camgrid, camhi)['SST'].isel(time=0)
   tszm = ts.zonal_mean(lat=sstlats)
   sstlat = tszm['latitudes'].isel(latitudes=tszm.argmax(dim='latitudes')).values
   sgnlat = np.sign(sstlat).astype(np.int_)
   print('%s: SST max lat %f has sign %d, SST value %.2f' % (sys.argv[0], sstlat, sgnlat, tszm.max(dim='latitudes').values))
   replace_or_add_variable(ptho, "sstmaxlat", sstlat)

if __name__ == '__main__':
   main()
