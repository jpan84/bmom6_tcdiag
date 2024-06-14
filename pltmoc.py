#plot merid overturning resid streamf based on https://mom6-analysiscookbook.readthedocs.io/en/latest/notebooks/Overturning.html
#import warnings
#warnings.filterwarnings('ignore')

import os
import numpy as np
import xarray as xr # requires >= 0.15.1
import pandas as pd
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt
from xoverturning import calcmoc

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0606dlyout/'
HISTS = 'ocn/hist/'
H0 = r'*mom6.hm_[0-9]*.nc'

def main():
   print('Opening history files...')
   ds = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H0)).rename(dict(zl='z_l', zi='z_i'))
   dsgrid = xr.open_dataset('/glade/derecho/scratch/jpan/b.e23.BMOM.f09_sx0.66av1.aqua.production.0515ctrl/run/INPUT/sx0.66av1_grid.nc')
   geo = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, '*static*.nc'))

   print(ds.dims)
   print(dsgrid.dims)

   moc = calcmoc(xr.merge([ds, geo]))#, dsgrid=dsgrid)
   print(moc)

if __name__ == '__main__':
   main()
