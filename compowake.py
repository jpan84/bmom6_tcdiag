### Take the strongest storms and plot their composite cold wakes. Relative to -7 to -3 day avg in a
### 10° box surrounding the point of max intensity
import sys
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import cftime
import datetime as dt
import matplotlib.pyplot as plt

### tempest traj output params
TRAJFILE = sys.argv[1] #output csv/parquet from traj_stats.py
BASEYR = 584 #year from which yeardelta is determined to stay within date bounds of pandas date objects

### hist file params
ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0815/'
HISTS = 'ocn/hist/'
H1 = '*mom6.hmd*'
DIRI = os.path.join(ARCHV, CASE, HISTS)
H1OFFSET_OCN = None #e.g., 6 if filename date is 01-01 and first timestamp is 06:00
STRF_OCN = lambda dtobj: f'*hmd_{dtobj.year:04}_{dtobj.month:02}_{dtobj.day:02}.nc'

### wake computation params
NTOP = 20 #number of strongest storms
LONBNDS = (-5, 5)
LATBNDS = (-5, 5)
AVBNDS = (-dt.timedelta(days=7), -dt.timedelta(days=3))
TBNDS = (-dt.timedelta(days=10), dt.timedelta(days=10))

omlvar = 'oml'
sstvar = 'tos'
budvars = {'hfsso': 'green', 'hflso': 'blue', 'rlntds': 'dimgray', 'rsntds': 'orange', 'Tadvconv': 'peru', 'Tdifconv': 'coral'}
budlabs = {'hfsso': 'SH', 'hflso': 'LH', 'rlntds': 'LW', 'rsntds': 'SW', 'Tadvconv': 'advection', 'Tdifconv': 'diffusion'}

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

   taxis = None
   omls = []
   ssts = []
   for ii, stm in enumerate(topstms.iterrows()):
      stm = stm[1]
      truedt = cftime.DatetimeNoLeap(stm['trueyr'], stm['month'], stm['day'], hour=stm['hour'])
      print(ii, truedt)

      filebnds = tuple([truedt + tbnd - dt.timedelta(hours=H1OFFSET_OCN) for tbnd in TBNDS])
      filebnds = tuple([glob.glob(os.path.join(DIRI, STRF_OCN(bnd)))[0] for bnd in filebnds]) #TODO: enforce that this be in bounds of the simulation dates
      toopen = h1s[(h1s >= filebnds[0]) & (h1s <= filebnds[1])]
      lonbnds = (stm['lon'] + LONBNDS[0], stm['lon'] + LONBNDS[1])
      latbnds = (stm['lat'] + LATBNDS[0], stm['lat'] + LATBNDS[1])
      ds = xr.open_mfdataset(toopen).sel(xh=slice(*lonbnds), yh=slice(*latbnds)) #TODO: use geolon/geolat to be more technically accurate

      omlser = ds[omlvar].mean(dim=['xh', 'yh'])
      sstser = ds[sstvar].mean(dim=['xh', 'yh'])
      omlref = omlser.sel(time=slice(truedt + AVBNDS[0], truedt + AVBNDS[1])).mean(dim='time')
      sstref = sstser.sel(time=slice(truedt + AVBNDS[0], truedt + AVBNDS[1])).mean(dim='time')

      budser = dict()
      for bk in budvars:
         budser[bk] = ds[bk].mean(dim=['xh', 'yh'])

      if not ii:
         taxis = [(tt - truedt).days + (tt - truedt).seconds/86400 for tt in omlser.time.values]
      #omls.append(omlser - omlref) #absolute anomaly
      omls.append(omlser / omlref * 100) #percentage
      ssts.append(sstser - sstref)

      #plt.plot(taxis, sstser)
      #plt.show()

   filoargs = (TRAJFILE.split('.')[0], NTOP, LONBNDS[1] - LONBNDS[0], LATBNDS[1] - LATBNDS[0])

   sstcomp = np.array([da.values for da in ssts]).mean(axis=0)
   plt.plot(taxis, sstcomp, color='red')
   plt.axvline(x=0, linestyle='--', color='black', linewidth=0.7)
   plt.xlabel('Day relative to max strength')
   plt.ylabel('SST relative to days %d to %d [K]' % (AVBNDS[0].days, AVBNDS[1].days))
   plt.title('Composite SST for top %d storms, lat %s, lon %s' % (NTOP, str(LATBNDS), str(LONBNDS)))
   ax2 = plt.gca().twinx()
   budlines = dict()
   for bk in budvars:
      budlines[bk] = ax2.plot(taxis, budser[bk].values, color=budvars[bk], label=budlabs[bk], linewidth=1)
   ###remove select lines
   for bk in ['Tadvconv', 'hfsso', 'rlntds']:
      budlines[bk][0].remove()
   ###
   ax2.legend()
   ax2.set_ylabel('Energy budget [W m-2]')
   plt.savefig('%s_%d_%dx%d_SSTwake.png' % filoargs, bbox_inches='tight')
   plt.close()

   omlcomp = np.array([da.values for da in omls]).mean(axis=0)
   plt.plot(taxis, omlcomp)
   plt.axvline(x=0, linestyle='--', color='black', linewidth=0.7)
   plt.xlabel('Day relative to max strength')
   plt.ylabel('OML relative to days %d to %d [%%]' % (AVBNDS[0].days, AVBNDS[1].days))
   plt.title('Composite OML for top %d storms, lat %s, lon %s' % (NTOP, str(LATBNDS), str(LONBNDS)))
   plt.savefig('%s_%d_%dx%d_%swake.png' % (*filoargs, omlvar))
   plt.close()


if __name__ == '__main__':
   main()
