#Joshua Pan Feb 2024
#Compute stats on the trajectories.txt output by par_track_aqua.sh (TempestExtremes)

#TODO: genesis lat
#TODO: termination lat
#TODO: check zonal sym
#TODO: peak winds
#TODO: cats
#TODO: SLP

import sys
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

MS2KT = 1.9438

SSHS = np.array([0, 17, 33, 43, 49, 58, 70, 9999])
COLS = [None, 'lonint', 'latp90', 'lon', 'lat', 'pres', 'wspd', 'year', 'month', 'day', 'hour']
TYPS = [None, int, int, float, float, float, float, int, int, int, int]
CTYP = {COLS[i+1]: TYPS[i+1] for i in range(len(COLS[1:]))}

XLIMS = dict()
YLIMS = dict()

def main(FN):
   print('Parsing files...')
   print(FN)
   f = open(FN, 'r')
   lns = f.readlines()
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
   #ids = list(range(len(trajs)))
   #print(ids)

   dfs = []
   for tn, tr in enumerate(trajs):
      #print(tr)
      stat = {var: [] for var in COLS[1:]}
      for ti in tr: #loop thru timesteps in trajectory
         [stat[vv].append(ti[jj + 1]) for jj, vv in enumerate(COLS[1:])]
      tcdf = pd.DataFrame.from_dict(stat).astype(CTYP)
      tcdf = tcdf.assign(mdstr=tcdf[COLS[-3:-1]].astype(str).agg('-'.join, axis=1))
      #dfs.append(tcdf.assign(dt=tcdf['dtstr'].astype(dt.datetime)))
      tcdf = tcdf.assign(doy=tcdf['mdstr'].apply(lambda x: int(dt.datetime.strptime(x, '%m-%d').strftime('%j'))))
      tcdf['yeardelta'] = tcdf['year'] % 584
      tcdf['year'] = tcdf['yeardelta'] + 1678
      tcdf = tcdf.assign(dt=pd.to_datetime(tcdf[COLS[-4:]]))
      #print(tcdf)
      print(tn)
      tcdf['stmnum'] = tn
      dfs.append(tcdf)

   print('Computing stats...')
   print(dfs[0])
   tc_stats = dict()
   tc_stats['genday'] = [tc['doy'].iloc[0] for tc in dfs]
   tc_stats['dura'] = [((tc['dt'].iloc[-1] - tc['dt'].iloc[0]) + (tc['dt'].iloc[1] - tc['dt'].iloc[0])).days for tc in dfs] #last-1st time plus 1 timestep
   tc_stats['pmins'] = [tc['pres'].min() for tc in dfs]
   tc_stats['genlat'] = [tc['lat'].iloc[0] for tc in dfs]
   tc_stats['lyslat'] = [tc['lat'].iloc[-1] for tc in dfs]
   tc_stats['maxu'] = [tc['wspd'].max() for tc in dfs]
   tc_stats['ace'] = [1e-4 * ((MS2KT * tc['wspd'])**2).sum() for tc in dfs]
   tc_stats['genlon'] = [tc['lon'].iloc[0] for tc in dfs]
   tc_stats['lyslon'] = [tc['lon'].iloc[-1] for tc in dfs]

   #cliplat = lambda lats: [np.nan if l < -15 else l for l in lats]
   #tc_stats['genlat'] = cliplat(tc_stats['genlat'])
   #tc_stats['lyslat'] = cliplat(tc_stats['lyslat'])

   outdf = pd.concat(dfs)
   outdf.to_parquet('%s.parquet' % FN.split('-')[-1])
   outdf.to_csv('%s.csv' % FN.split('-')[-1])

   #FN = 'tcstats-QPC4ctrl'
   print('Hists and scatters...')
   for k in tc_stats:
      print(k)
      plt.hist(tc_stats[k], density=True, bins=15, edgecolor='black')
      plt.title(k)
      plt.ylabel('Probability density')
      if k in XLIMS:
         plt.xlim(*XLIMS[k])
         plt.ylim(*YLIMS[k])
      plt.savefig('%s/%s.png' % (FN.split('-')[-1], k))
      plt.close()
      for depvar in tc_stats:
         plt.scatter(tc_stats[k], tc_stats[depvar])
         plt.xlabel(k)
         plt.ylabel(depvar)
         plt.savefig('%s/scat_%s_%s.png' % (FN.split('-')[-1], k, depvar))
         plt.close()
   
   print('%s done.' % sys.argv[0])

def mps2cat(wspd):
   dif = wspd - SSHS
   return np.argmax(dif < 0) - 2

if __name__ == '__main__':
   FN = sys.argv[1]
   main(FN)
