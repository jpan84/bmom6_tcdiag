#Joshua Pan Mar 2025
#Tweaked to use unstructured grid
#Compute stats on the trajectories.txt output by par_track_aqua.ne120.CZ.sh (TempestExtremes)

#TODO: genesis lat
#TODO: termination lat
#TODO: check zonal sym
#TODO: peak winds
#TODO: cats
#TODO: SLP

import os
import sys
import pickle
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
from scipy import stats

MS2KT = 1.9438

SSHS = np.array([0, 17, 33, 43, 49, 58, 70, 9999])
COLS = [None, 'ncol', 'lon', 'lat', 'pres', 'wspd', 'year', 'month', 'day', 'hour']
TYPS = [None, int, float, float, float, float, int, int, int, int]
CTYP = {COLS[i+1]: TYPS[i+1] for i in range(len(COLS[1:]))}

XLIMS = dict(pmins=(8.5e4,1.01e5), ace=(0,80), maxu=(15,120), genlon=(0,360), genlat=(-40, 40))
YLIMS = dict(pmins=(0, 6e-2))
clabelkwargs = {'inline': 1, 'fontsize': 10, 'colors': 'black', 'fmt': '%.1e'}

def main(FN, CTRL=None):
   DOUT = FN.split('.')[-1]
   if not os.path.exists(DOUT):
      os.makedirs(DOUT)

   print('Parsing files...')
   print(FN)
   f = open(FN, 'r')
   lns = f.readlines()
   tbl = [l.replace('\n', '').split('\t') for l in lns]

   ctrl_kde = None
   if not CTRL is None:
      with open(CTRL, 'rb') as handle:
         ctrl_kde = pickle.load(handle)

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
      #dfs.append(tcdf.assign(dt=tcdf['dtstr'].astype(dt.datetime)))
      
      tcdf = tcdf.assign(mdstr=tcdf[COLS[-3:-1]].astype(str).agg('-'.join, axis=1))
      tcdf = tcdf.assign(doy=tcdf['mdstr'].apply(lambda x: int(dt.datetime.strptime(x, '%m-%d').strftime('%j'))))
      '''
      tcdf['trueyr'] = tcdf['year']
      tcdf['yeardelta'] = tcdf['year'] % 584
      tcdf['year'] = tcdf['yeardelta'] + 1678
      tcdf = tcdf.assign(dt=pd.to_datetime(tcdf[COLS[-4:]]))
      '''
      #turn date cols into zero-padded strings
      tcdf = tcdf.assign(year=tcdf['year'] + 2000)#.astype(str).str.zfill(4))
      tcdf = tcdf.assign(month=tcdf['month'].astype(str).str.zfill(2))
      tcdf = tcdf.assign(day=tcdf['day'].astype(str).str.zfill(2)) 
      #print(tcdf[COLS[-4:]])
      tcdf = tcdf.assign(dt=pd.to_datetime(tcdf[COLS[-4:-1]]) + pd.to_timedelta(tcdf['hour'], unit='h'))

      #resample winds to 6H for ACE
      tcdf = tcdf.set_index('dt')
      tcdf = tcdf.resample('6H').agg({col: 'first' for col in tcdf.columns})
      tcdf['wspd'] = tcdf['wspd'].interpolate(method='linear')
      
      #print(tcdf)
      print(tn)
      tcdf['stmnum'] = tn
      dfs.append(tcdf)

   print('Computing stats...')
   print(dfs[:1])
   tc_stats = dict()
   tc_stats['genday'] = [tc['doy'].iloc[0] for tc in dfs]
   #tc_stats['dura'] = [((tc['dt'].iloc[-1] - tc['dt'].iloc[0]) + (tc['dt'].iloc[1] - tc['dt'].iloc[0])).days for tc in dfs] #last-1st time plus 1 timestep
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
   outdf.to_parquet('%s.parquet' % DOUT)
   outdf.to_csv('%s.csv' % DOUT)

   #FN = 'tcstats-QPC4ctrl'
   print('Hists and scatters...')
   plt.rc('font', size=20)
   kde_dict = dict()
   for k in tc_stats:
      print(k)
      if k == 'pmins':
         plt.hist([pval for pval in tc_stats[k]], density=True, bins=15, edgecolor='black') #pval / 100 for hPa
      else:
         plt.hist(tc_stats[k], density=True, bins=15, edgecolor='black')
      if k == 'ace':
         plt.title('%s (total=%.1f)' % (k, np.sum(tc_stats[k])))
      else:
         plt.title(k)
      plt.ylabel('Probability density')
      if k in XLIMS:
         plt.xlim(*XLIMS[k])
      if k in YLIMS:
         plt.ylim(*YLIMS[k])
      plt.savefig('%s/%s.png' % (DOUT, k), bbox_inches='tight')
      plt.close()
      for depvar in tc_stats:
         if depvar == k:
            continue
         # Start with a square Figure.
         fig = plt.figure(figsize=(6, 6))
         # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
         # the size of the marginal Axes and the main Axes in both directions.
         # Also adjust the subplot parameters for a square plot.
         gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                               left=0.1, right=0.9, bottom=0.1, top=0.9,
                               wspace=0.05, hspace=0.05)
         # Create the Axes.
         ax = fig.add_subplot(gs[1, 0])
         ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
         ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
         # Draw the scatter plot and marginals.
         kde_k_dep = scatter_hist(tc_stats[k], tc_stats[depvar], ax, ax_histx, ax_histy, k, depvar, ctrl_kde=ctrl_kde)
         kde_dict[(k, depvar)] = kde_k_dep
         ax.set_xlabel(k)
         ax.set_ylabel(depvar)
         #plt.scatter(tc_stats[k], tc_stats[depvar])
         #plt.xlabel(k)
         #plt.ylabel(depvar)
         plt.savefig('%s/scat_%s_%s.png' % (DOUT, k, depvar), bbox_inches='tight')
         plt.close()

   with open('%s.kde.pickle' % DOUT, 'wb') as handle:
      pickle.dump(kde_dict, handle, protocol=-1)
   print('%s done.' % sys.argv[0])

def scatter_hist(x, y, ax, ax_histx, ax_histy, xname, yname, ctrl_kde=None):
   # no labels
   ax_histx.tick_params(axis="x", labelbottom=False, labelleft=False)
   ax_histy.tick_params(axis="y", labelleft=False, labelbottom=False)

   #kde beneath the scatter
   points = np.vstack([x, y])
   kde = stats.gaussian_kde(points)
   xmin, xmax, ymin, ymax = min(x), max(x), min(y), max(y)
   if xname in XLIMS:
      xmin, xmax = XLIMS[xname]
   if yname in XLIMS:
      ymin, ymax = XLIMS[yname]
   xgrid, ygrid = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
   pos = np.vstack([xgrid.ravel(), ygrid.ravel()])
   z = kde(pos).reshape(xgrid.shape)

   if ctrl_kde is None:
      csf = ax.contourf(xgrid, ygrid, z)
      cs = ax.contour(xgrid, ygrid, z, colors='black')
      ax.clabel(cs, **clabelkwargs)
   else:
      zctrl = ctrl_kde[(xname, yname)](pos).reshape(xgrid.shape)
      csf = ax.contourf(xgrid, ygrid, z - zctrl, cmap='bwr')
      cs = ax.contour(xgrid, ygrid, zctrl, colors='black')
      plt.colorbar(csf)

   # the scatter plot:
   ax.scatter(x, y)

   ax_histx.hist(x, bins=20, edgecolor='black')
   ax_histy.hist(y, bins=20, edgecolor='black', orientation='horizontal')
   ax_histx.set_yticks([])
   ax_histy.set_xticks([])

   return kde

def mps2cat(wspd):
   dif = wspd - SSHS
   return np.argmax(dif < 0) - 2

if __name__ == '__main__':
   FN = sys.argv[1]
   CTRL = None #pickled dict of kde objects
   if len(sys.argv) > 2:
      CTRL = sys.argv[2]
   main(FN, CTRL=CTRL)
