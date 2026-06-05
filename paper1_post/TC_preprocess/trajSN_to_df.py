#Joshua Pan
#parse TempestExtremes StitchNode files into pandas dataframe/csv/parquet

import os
import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
from paths import ARCHRT, CASENAMES
import numpy as np
import pandas as pd
import datetime as dt

TRAJFILES = ['trajectories.txt.' + lc for lc in CASENAMES]

SSHS = np.array([0, 17, 33, 43, 49, 58, 70, 9999])
COLS = [None, 'ncol', 'lon', 'lat', 'pres', 'wspd', 'year', 'month', 'day', 'hour']
TYPS = [None, int, float, float, float, float, int, int, int, int]
CTYP = {COLS[i+1]: TYPS[i+1] for i in range(len(COLS[1:]))}

def main(FN, CTRL=None):
   print('Parsing files...')
   print(FN)
   f = open(FN, 'r')
   lns = f.readlines()
   if len(lns) == 0:
      return
   tbl = [l.replace('\n', '').split('\t') for l in lns]

   trajs = []
   curr_ts = []
   for ln in tbl:
      if ln[0] == 'start':
         if len(curr_ts):
            trajs.append(curr_ts)
            curr_ts = []
      else:
         curr_ts.append(ln)
   trajs.append(curr_ts)

   print(len(trajs))

   dfs = []
   for tn, tr in enumerate(trajs):
      stat = {var: [] for var in COLS[1:]}
      for ti in tr: #loop thru timesteps in trajectory
         [stat[vv].append(ti[jj + 1]) for jj, vv in enumerate(COLS[1:])]
      tcdf = pd.DataFrame.from_dict(stat).astype(CTYP)
      
      tcdf = tcdf.assign(mdstr=tcdf['month'].astype(str) + '-' + tcdf['day'].astype(str))
      #tcdf = tcdf.assign(mdstr=tcdf[COLS[-3:-1]].astype(str).agg('-'.join, axis=1))
      tcdf = tcdf.assign(doy=tcdf['mdstr'].apply(lambda x: int(dt.datetime.strptime(x, '%m-%d').strftime('%j'))))

      #turn date cols into zero-padded strings
      tcdf = tcdf.assign(year=tcdf['year'] + 2000)#.astype(str).str.zfill(4))
      tcdf = tcdf.assign(month=tcdf['month'].astype(str).str.zfill(2))
      tcdf = tcdf.assign(day=tcdf['day'].astype(str).str.zfill(2)) 
      tcdf = tcdf.assign(dt=pd.to_datetime(tcdf[COLS[-4:-1]]) + pd.to_timedelta(tcdf['hour'], unit='h'))

      tcdf = tcdf.set_index('dt')
      tcdf = tcdf.resample('6h').agg({col: 'first' for col in tcdf.columns})

      print(tn)
      tcdf['stmnum'] = tn
      tcdf['isgen'] = (tcdf.index == tcdf.index[0])
      tcdf['islys'] = (tcdf.index == tcdf.index[-1])
      dfs.append(tcdf)

   outdf = pd.concat(dfs)
   outdf.to_parquet('%s.parquet' % FN)
   outdf.to_csv('%s.csv' % FN)

if __name__ == '__main__':
   for tfl in TRAJFILES:
      main(tfl)
