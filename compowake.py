### Take the strongest storms and plot their composite cold wakes. Relative to -7 to -3 day avg in a
### 10° box surrounding the point of max intensity
import os
import numpy as np
import xarray as xr
import pandas as pd
import cftime
import datetime as dt
import matplotlib.pyplot as plt

### tempest traj output params
TRAJFILE = 'set4.parquet' #output csv/parquet from traj_stats.py

### hist file params
ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0711dlyout/'
HISTS = 'ocn/hist/'
H1OFFSET = None

### wake computation params
NTOP = 5 #number of strongest storms
LONBNDS = (-10, 10)
LATBNDS = (-10, 10)
AVBNDS = (-dt.timedelta(days=7), -dt.timedelta(days=3))

omlvar = 'oml'
sstvar = 'tos'

def main():
   tcdf = pd.read_parquet(TRAJFILE)
   maxu = tcdf['wspd'].groupby(tcdf['stmnum']).max()
   print(maxu)

   peaks = tcdf.groupby('stmnum').apply(lambda x: x.loc[x['wspd'].idxmax()])
   print(peaks)

   exit()

   print('Opening history files...')
   ds = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H1)).sel(time=slice(*TBNDS), xh=slice(*LONBNDS), yh=slice(*LATBNDS))

   trc = ds[var].mean(dim=['xh', 'yh'])
   taxis = [(tt - pltdate).days + (tt - pltdate).seconds/86400 for tt in trc.time.values]

   print(taxis)

   plt.plot(taxis, trc.values)
   plt.title('%s %s' % (str(pltdate), var))
   plt.show()

if __name__ == '__main__':
   main()
