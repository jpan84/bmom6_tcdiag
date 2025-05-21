#construct a pandas dataframe of un/seed events
import os
import glob
import re
import numpy as np
import pandas as pd
from datetime import datetime as dt
'''
import uxarray as ux

import holoviews as hv
from holoviews.operation import contours as hvcontours
from bokeh.models import FixedTicker
import geoviews.feature as gf
import cartopy.crs as ccrs
'''

alias = '250415_unseed_production'
dfcols = ['dt', 'sstlat', 'sstval', 'tid', 'clat', 'clon', 'psmin', 'rp', 'dp', 'zp', 'exppr']
logpaths = sorted(glob.glob('/glade/u/home/jpan/work/MOM6_CASEDIRS/%s.out*' % alias))

'''
BBOX_DEG = 10.
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250321_unseed_400.15.8e-5'
OUTDIR = 'bef_aft.250321_unseed_400.15.8e-5'
seedlog = open('/glade/u/home/jpan/work/MOM6_CASEDIRS/250321_unseed_400.15.8e-5.out1', 'r')
RESTDIR = '/glade/derecho/scratch/jpan/archive/%s/rest/%s/'
HTRACKSTR = '*.cam.h1i.*.nc'
grid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
'''

ps_pattern_1000 = r"^\[1\d{5}\."
ps_pattern_900 = r"^9\d{4}\."

'''
dpmin = 0 #3 #hPa
dTmax = 0 #-1 #K
CEN_DEG = 3 #find largest dp and dT within this many GCD of clon,clat

Uclim = 20

lev = None
'''

def main():
   #print('\n\n\nStarting with dp, dT thresholds\n', dpmin, dTmax)
   #global lev
   #if not os.path.exists(OUTDIR):
   #   os.makedirs(OUTDIR)

   #print(list(open(logpaths[0], 'r')))
   seedlog = list(open(logpaths[0], 'r'))
   for pt in logpaths[1:]:
      seedlog.extend(list(open(pt, 'r')))
   #print(seedlog)
   #exit()

   outdf = pd.DataFrame(columns=dfcols)
   mkplt = False
   clat, clon, psmin, tid, sstlat, sstval, dtstr = tuple([np.nan for _ in range(7)])
   settings = []
   #dtstr, orinc, modnc = None, None, None
   prevln = ''
   for ln in seedlog:
      spl = ln.strip('\n').split(' ')
      if not len(spl[0]):
         continue
      if 'not in same hemisphere' in ln and ln != prevln:
         psmin, settings = np.nan, []
         mkplt = True
      if spl[0] == 'find-sst-max.py:':
         sstlat = float(spl[4])
         sstval = float(spl[-1])
      if spl[0] == 'sedding':
         #211 sedding lat/lon: -8.347259 48.750000
         clon = float(spl[-1])
         clat = float(spl[-2])
      elif spl[:2] == ['Debugging', 'restart']:
         #212 Debugging restart file /glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc at time 5407747
         tid = int(spl[-1])
      elif spl[0][0] == "'":
         #200 '/glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc' -> '/glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc.ORIG.nc'
         dtstr = re.search(r"(\d{4}-\d{2}-\d{2}-\d{5})", spl[0]).group()
         print(dtstr)
         #modnc = os.path.basename(spl[0].strip("'"))
         #orinc = os.path.basename(spl[-1].strip("'"))
         #print(modnc, orinc)
      elif spl[0][0] == '[':
         if re.match(ps_pattern_1000, spl[0]): 
            psmin = float(spl[0].strip('[')) / 100
            #mkplt = True #make plot now that we have all the info for 1 storm
         if spl[0] == '[' and re.match(ps_pattern_900, spl[1]):
            psmin = float(spl[1]) / 100
            #mkplt = True
      if 'Final BEST SETTINGS' in ln:
         settings = [float(num) for num in re.compile(r'\d+\.\d+').findall(ln)]
         mkplt = True
      if mkplt:
         newrow = [dt.strptime('-'.join(dtstr.split('-')[:-1]), '%Y-%m-%d'), sstlat, sstval, tid, clat, clon, psmin, *settings]
         ser = pd.Series(newrow, index=outdf.columns[:len(newrow)])
         outdf.loc[len(outdf)] = ser
         mkplt = False
      prevln = ln

   outdf.to_parquet('seed_stats/%s_unseed_events.parquet' % alias)
   outdf.to_csv('seed_stats/%s_unseed_events.csv' % alias)

if __name__ == '__main__':
   main()
