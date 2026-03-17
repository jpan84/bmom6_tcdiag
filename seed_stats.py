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

alias = '251229_seed_match'
dfcols = ['dt', 'sstlat', 'sstval', 'tid', 'clat', 'clon', 'psamb', 'dp [hPa]', 'rmw_target', 'rmw_final', 'rp']
logpaths = sorted(glob.glob('/glade/u/home/jpan/work/MOM6_CASEDIRS/%s.out*' % alias))

def main():
   seedlog = list(open(logpaths[0], 'r'))
   for pt in logpaths[1:]:
      seedlog.extend(list(open(pt, 'r')))

   outdf = pd.DataFrame(columns=dfcols)
   mkplt = False
   clat, clon, psamb, tid, sstlat, sstval, dtstr, dp, rmwt, rmwf, rp = tuple([np.nan for _ in range(11)])
   prevln = ''
   for ln in seedlog:
      spl = ln.strip('\n').split(' ')
      if not len(spl[0]):
         continue

      if spl[0] == 'find-sst-max.py:':
         sstlat = float(spl[4])
         sstval = float(spl[-1])
      elif spl[0][0] == "'":
         #200 '/glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc' -> '/glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc.ORIG.nc'
         dtstr = re.search(r"(\d{4}-\d{2}-\d{2}-\d{5})", spl[0]).group()
         print(dtstr)
      elif spl[:2] == ['Debugging', 'restart']:
         #212 Debugging restart file /glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc at time 5407747
         tid = int(spl[-1])
      elif spl[0] == 'lat/lon:':
         clat, clon = float(spl[-2]), float(spl[-1])
      elif spl[0] == 'ambient':
         psamb = float(spl[-1])
      elif spl[0] == 'iter:':
         #print(spl)
         rmwt = float(spl[-4])
      elif spl[:2] == ['random', 'minp:']:
         dp = float(spl[-1])
      elif spl[0] == 'rp:':
         rp = float(spl[1])
         mkplt = True
      elif prevln[:5] == 'iter:' and spl[0] != 'iter:':
         rmwf = float(prevln.strip('\n').split(' ')[7])

      if mkplt:
         newrow = [dt.strptime('-'.join(dtstr.split('-')[:-1]), '%Y-%m-%d'), sstlat, sstval, tid, clat, clon, psamb, dp, rmwt, rmwf, rp]
         #print(outdf.columns)
         #print(newrow)
         ser = pd.Series(newrow, index=outdf.columns[:len(newrow)])
         outdf.loc[len(outdf)] = ser
         mkplt = False

      prevln = ln

   outdf.to_parquet('seed_stats/%s_seed_events.parquet' % alias)
   outdf.to_csv('seed_stats/%s_seed_events.csv' % alias)

if __name__ == '__main__':
   main()
