#plot oht and explore variability
#import warnings
#warnings.filterwarnings('ignore')

import os
import numpy as np
import xarray as xr # requires >= 0.15.1
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0530ywbranch/'
HISTS = 'ocn/hist/'
H0 = r'*mom6.hm_[0-9]*.nc'

def main():
   print('Opening history files...')
   ds = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H0)).mean(dim='time').sum(dim='xh')

   oht = ds.T_ady_2d + ds.T_diffy_2d
   plt.plot(oht.yq, oht.values)
   plt.show()

if __name__ == '__main__':
   main()
