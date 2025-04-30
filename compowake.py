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
TRAJFILE = sys.argv[1] #output parquet from traj_stats.py

### hist file params
ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1'
HISTS = 'ocn/hist/'
H1 = '*mom6.sfc*'
STATIC = '/glade/derecho/scratch/jpan/%s/run/%s.mom6.static.nc' % CASE
#TODO: account for choice of monthly or daily sfc history file freq
DIRI = os.path.join(ARCHV, CASE, HISTS)
H1OFFSET_OCN = None #e.g., 6 if filename date is 01-01 and first timestamp is 06:00
STRF_OCN = lambda dtobj: f'*sfc.{dtobj.year:04}-{dtobj.month:02}-{dtobj.day:02}.nc'

### wake computation params
NTOP = 10 #number of strongest storms
LONBNDS = (-5, 5)
LATBNDS = (-5, 5)
AVBNDS = (-dt.timedelta(days=7), -dt.timedelta(days=3))
TBNDS = (-dt.timedelta(days=10), dt.timedelta(days=10))

omlvar = 'oml'
sstvar = 'tos'
budvars = {'hfsso': 'green', 'hflso': 'blue', 'rlntds': 'dimgray', 'rsntds': 'orange', 'Tadvconv': 'peru', 'Tdifconv': 'coral', 'hfds': 'black'}
budlabs = {'hfsso': 'SH', 'hflso': 'LH', 'rlntds': 'LW', 'rsntds': 'SW', 'Tadvconv': 'advection', 'Tdifconv': 'diffusion', 'hfds': 'sfctot'}

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
   #print(type(topstms.iloc[0]['dt']))
   geogrid = xr.open_dataset(STATIC)

   print('Computing composites...')
   taxis = None
   omls = []
   ssts = []
   for ii, stm in enumerate(topstms.iterrows()):
      stm = stm[1]
      truedt = cftime.DatetimeNoLeap(int(stm['year']) - 2000, int(stm['month']), int(stm['day']), hour=int(stm['hour']))
      print(ii, truedt)

      lonbnds = ((stm['lon'] + LONBNDS[0]) % 360., (stm['lon'] + LONBNDS[1]) % 360.)
      latbnds = (stm['lat'] + LATBNDS[0], stm['lat'] + LATBNDS[1])
      #TODO: make wherelat/lon into a function to account for q points too
      wherelat = (geogrid.geolat >= latbnds[0]) & (geogrid.geolat <= latbnds[1])
      wherelon = (geogrid.geolon >= lonbnds[0]) & (geogrid.geolon <= lonbnds[1])
      if lonbnds[1] < lonbnds[0]:
         wherelon = (geogrid.geolon >= lonbnds[0]) | (geogrid.geolon <= lonbnds[1])
      selarea = lambda da: da.where(wherelat & wherelon, drop=True)
      sellat = selarea(geogrid.geolat)

      #open files with the desired times
      filebnds = tuple([truedt + tbnd - dt.timedelta(hours=H1OFFSET_OCN) for tbnd in TBNDS])
      #print(STRF_OCN(filebnds[0]))
      filebnds = tuple([glob.glob(os.path.join(DIRI, STRF_OCN(bnd)))[0] for bnd in filebnds]) #TODO: enforce that this be in bounds of the simulation dates
      toopen = h1s[(h1s >= filebnds[0]) & (h1s <= filebnds[1])]
      ds = xr.open_mfdataset(toopen)#.sel(xh=slice(*lonbnds), yh=slice(*latbnds)) #TODO: use geolon/geolat to be more technically accurate

      #print(ds[sstvar].where(wherelat & wherelon, drop=True))

      omlser = latlon_avg(selarea(ds[omlvar]), sellat)
      sstser = latlon_avg(selarea(ds[sstvar]), sellat)
      #print(sstser)
      omlref = omlser.sel(time=slice(truedt + AVBNDS[0], truedt + AVBNDS[1])).mean(dim='time') #reference values -7 to -3 days
      sstref = sstser.sel(time=slice(truedt + AVBNDS[0], truedt + AVBNDS[1])).mean(dim='time')

      budser = dict()
      for bk in budvars:
         budser[bk] = latlon_avg(selarea(ds[bk]), sellat) #TODO: area weight the mean

      if not ii:
         taxis = [(tt - truedt).days + (tt - truedt).seconds/86400 for tt in omlser.time.values]
      #omls.append(omlser - omlref) #absolute anomaly
      omls.append(omlser / omlref * 100) #percentage
      ssts.append(sstser - sstref)
      #print(ssts)

      #plt.plot(taxis, sstser)
      #plt.show()

   filoargs = (os.path.splitext(os.path.basename(TRAJFILE))[0], NTOP, LONBNDS[1] - LONBNDS[0], LATBNDS[1] - LATBNDS[0])

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
   ###remove select lines for presentation slides
   #for bk in ['hfsso', 'rlntds', 'hflso', 'rsntds', 'Tadvconv', 'hfds']:
   #   budlines[bk][0].remove()
   ###
   #mybud = sum([budser[bk] for bk in budser if bk != 'hfds'])
   #print(mybud)
   #ax2.plot(taxis, mybud - budser['hfds'], color = 'purple', label='resid', linewidth=1.5) #this residual ended up closely matching Tadvconv
   ax2.plot(taxis, budser['Tadvconv'] + budser['hfds'], color = 'purple', label='adv+sfc')
   ax2.legend()
   ax2.set_ylabel('Energy budget [W m-2]')
   plt.savefig('%s_%d_%dx%d_SSTwake.png' % filoargs, bbox_inches='tight')
   plt.close()

   #TODO: plot something with wind stress
   omlcomp = np.array([da.values for da in omls]).mean(axis=0)
   plt.plot(taxis, omlcomp)
   plt.axvline(x=0, linestyle='--', color='black', linewidth=0.7)
   plt.xlabel('Day relative to max strength')
   plt.ylabel('OML relative to days %d to %d [%%]' % (AVBNDS[0].days, AVBNDS[1].days))
   plt.title('Composite OML for top %d storms, lat %s, lon %s' % (NTOP, str(LATBNDS), str(LONBNDS)))
   plt.savefig('%s_%d_%dx%d_%swake.png' % (*filoargs, omlvar))
   plt.close()

def latlon_avg(da, latcoord, latdim='yh', londim='xh'):
   wgts = np.cos(np.deg2rad(latcoord))
   return ((da * wgts).sum(dim=latdim) / wgts.sum(dim=latdim)).mean(dim=londim)

if __name__ == '__main__':
   main()
