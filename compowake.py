### Take the strongest storms and plot their composite cold wakes. Relative to -7 to -3 day avg in a
### 10° box surrounding the point of max intensity
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import cftime
import datetime as dt
import matplotlib.pyplot as plt

### tempest traj output params
TRAJFILE = 'set4.parquet' #output csv/parquet from traj_stats.py
BASEYR = 584 #year from which yeardelta is determined to stay within date bounds of pandas date objects

### hist file params
ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0711dlyout/'
HISTS = 'ocn/hist/'
H1 = '*mom6.hmd*'
DIRI = os.path.join(ARCHV, CASE, HISTS)
H1OFFSET_OCN = None #e.g., 6 if filename date is 01-01 and first timestamp is 06:00
STRF_OCN = lambda dtobj: f'*hmd_{dtobj.year:04}_{dtobj.month:02}_{dtobj.day:02}.nc'

### wake computation params
NTOP = 5 #number of strongest storms
LONBNDS = (-10, 10)
LATBNDS = (-10, 10)
AVBNDS = (-dt.timedelta(days=7), -dt.timedelta(days=3))
TBNDS = (-dt.timedelta(days=10), dt.timedelta(days=10))

omlvar = 'oml'
sstvar = 'tos'

def main():
   tcdf = pd.read_parquet(TRAJFILE)
   maxu = tcdf['wspd'].groupby(tcdf['stmnum']).max()
   print(maxu)

   peaks = tcdf.groupby('stmnum').apply(lambda x: x.loc[x['wspd'].idxmax()])
   print(peaks)

   topstms = peaks.nlargest(NTOP, 'wspd')
   print(topstms)

   print('Getting list of history files...')
   h1s = sorted(glob.glob(os.path.join(DIRI, H1)))
   h1s = np.array(h1s)
   rndds = xr.open_dataset(h1s[0])
   #print(h1s[0])
   H1OFFSET_OCN = rndds.time.values[0].hour
   print(type(topstms.iloc[0]['dt']))

   print('Computing composites...')
   '''
   #trueyr = BASEYR + topstms.iloc[0]['yeardelta']
   truedt = cftime.DatetimeNoLeap(topstms.iloc[0]['trueyr'], topstms.iloc[0]['month'], topstms.iloc[0]['day'], hour=topstms.iloc[0]['hour'])
   print(truedt)

   filebnds = tuple([truedt + tbnd - dt.timedelta(hours=H1OFFSET_OCN) for tbnd in TBNDS])
   filebnds = tuple([glob.glob(os.path.join(DIRI, STRF_OCN(bnd)))[0] for bnd in filebnds]) #TODO: enforce that this be in bounds of the simulation dates
   toopen = h1s[(h1s >= filebnds[0]) & (h1s <= filebnds[1])]
   ds = xr.open_mfdataset(toopen)
   print(toopen)
   '''

   for ii, stm in enumerate(topstms.iterrows()):
      stm = stm[1]
      truedt = cftime.DatetimeNoLeap(stm['trueyr'], stm['month'], stm['day'], hour=stm['hour'])
      print(ii, truedt)

      filebnds = tuple([truedt + tbnd - dt.timedelta(hours=H1OFFSET_OCN) for tbnd in TBNDS])
      filebnds = tuple([glob.glob(os.path.join(DIRI, STRF_OCN(bnd)))[0] for bnd in filebnds]) #TODO: enforce that this be in bounds of the simulation dates
      toopen = h1s[(h1s >= filebnds[0]) & (h1s <= filebnds[1])]
      ds = xr.open_mfdataset(toopen)

      print(ds)










   exit()
   truedt = lambda df: cftime.DatetimeNoLeap(df[['trueyr', 'month', 'day', 'hour']])
   print(truedt(topstms))
   exit()

   truedt = topstms.iloc[0]['dt']
   truedt = truedt + dt.timedelta(days=365*int(-truedt.year + BASEYR + topstms.iloc[0]['yeardelta']))
   print(truedt)
   exit()

   ds = xr.open_dataset(os.path.join(ARCHV, CASE, HISTS, H1)).sel(time=slice(*TBNDS), xh=slice(*LONBNDS), yh=slice(*LATBNDS))

   exit()
   trc = ds[var].mean(dim=['xh', 'yh'])
   taxis = [(tt - pltdate).days + (tt - pltdate).seconds/86400 for tt in trc.time.values]

   print(taxis)

   plt.plot(taxis, trc.values)
   plt.title('%s %s' % (str(pltdate), var))
   plt.show()

if __name__ == '__main__':
   main()
