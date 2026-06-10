#Joshua Pan Mar 2025
#Uses parquet output of traj_stats_bmom_SE.py

import os
import sys
import pickle
import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats

PLOTSEEDS = int(sys.argv[1])

MS2KT = 1.9438
a_e = 6.371e6 #m
LATLABS = np.arange(-90, 90.1, 30)
FIXYLIM = (0, 70)
DLAT = 1.5

SZNS = ['DJF', 'MAM', 'JJA', 'SON']
SZMOs = [{12, 1, 2}, {3, 4, 5}, {6, 7, 8}, {9, 10, 11}]

FILI = ['250702_unseed2hPa6m.csv.flagged.parquet', '250415_unseed.csv.flagged.parquet', '250417_ctrl.parquet', '251229_seedmatch.seedflagged.parquet', '250416_seed1x1.seedflagged.parquet'] #flagged by un/seed_success_rate.py
labels = ['unseed2', 'unseed', 'ctrl', 'mseed', 'seed']
events = [('250702_unseed_2hPa6m_unseed_events.parquet', 'us'), ('250415_unseed_production_unseed_events.parquet', 'us'), None, ('251229_seed_match_seed_events.parquet', 'sd'), ('250416_seed1x1_production_seed_events.parquet', 'sd')]
rename_dict = dict(clat='lat', clon='lon')
DOUT = '260415_density_5exp'

#FILI = ['250415_unseed_JJASOND.parquet', '250417_ctrl_JJASOND.parquet']
#labels = ['unseed', 'ctrl']
#DOUT = '250423_density_JJASOND'

nhmask, shmask = None, None #boolean mask is set in make_bin_grid_1d to help hemi selection

def main():
   #extidx = FILI.rfind('.')
   #DOUT = FILI[:extidx] + '_density'
   if not os.path.exists(DOUT):
      os.makedirs(DOUT)
   open(os.path.join(DOUT, 'sznl_climo.csv'), 'w').close() #reset the seasonal-total csv file

   dfs = [pd.read_parquet(f) for f in FILI]
   for df in dfs:
      df['max_lft_wspd'] = df.groupby('stmnum')['wspd'].transform('max')
      if type(df.index) != pd.core.indexes.datetimes.DatetimeIndex:
         df['dt'] = pd.to_datetime(df['dt'])
         df.set_index('dt', inplace=True)
   #TODO !1124: add a flag indicating the true genesis point before splitting into months
   print(dfs[0].index, type(dfs[0].index))
   #print(dfs[0].groupby(dfs[0].index.month))
   trange = [df.index.max() - df.index[0] for df in dfs]
   print(trange)
   totyrs = [tr.total_seconds() / 86400 / 365 for tr in trange]
   print(totyrs)
   #exit()

   evdfs = []
   if PLOTSEEDS:
      #usdf = pd.read_parquet(events[0]).rename(columns=rename_dict)
      #sddf = pd.read_parquet(events[-1]).rename(columns=rename_dict)
      evdfs = [pd.read_parquet(ev[0]).rename(columns=rename_dict) if ev is not None else None for ev in events]
      #print(usdf)
      #print(usdf['dt'])

   bininfo = make_bin_grid_1d(latbnds=(-90., 90.), dlat=DLAT)
   global nhmask, shmask
   nhmask = (bininfo[1] > 0) & (abs(bininfo[1]) > 1e-2) #omits bins centered at 0 deg lat
   shmask = (bininfo[1] < 0) & (abs(bininfo[1]) > 1e-2)
   print(nhmask, shmask)

   plt.rc('font', size=14)

   outdss = []
   for sznnm, sznmos in zip(SZNS, SZMOs):
      #print('\nWorking on season', sznnm)
      szn_dfs = [df[df.index.month.isin(sznmos)] for df in dfs]

      uniq = plot_lat_binned(szn_dfs, bininfo, totyrs, unique_tracks_of_storm_1d, 'unique storms', 'unique_track_dens_%.1f_%s.png' % (DLAT, sznnm))
      h6all = plot_lat_binned(szn_dfs, bininfo, totyrs, sixhrly_fixes_of_storm_1d, '6-hourly fixes', '6hr_track_dens_%.1f_%s.png' % (DLAT, sznnm))
      h6hurr = plot_lat_binned(szn_dfs, bininfo, totyrs, sixhrly_fixes_of_storm_1d, '6-hourly hurricane fixes', '6hr_track_dens_hurr_%.1f_%s.png' % (DLAT, sznnm), wspd_thresh=33)
      h6maj = plot_lat_binned(szn_dfs, bininfo, totyrs, sixhrly_fixes_of_storm_1d, '6-hourly major hurricane fixes', '6hr_track_dens_major_%.1f_%s.png' % (DLAT, sznnm), wspd_thresh=50)
      ace = plot_lat_binned(szn_dfs, bininfo, totyrs, bin_ace_of_storm_1d, 'ACE [$10^4$ kt$^2$] ', 'ace_binned_%.1f_%s.png' % (DLAT, sznnm))
      gen = plot_lat_binned(szn_dfs, bininfo, totyrs, genesis_pts_1d, 'genesis points (storm count)', 'gen_pt_dens_%.1f_%s.png' % (DLAT, sznnm))
      sd_gen = plot_lat_binned(szn_dfs, bininfo, totyrs, genesis_pts_1d, 'seeded genesis points', 'sd_gen_pt_dens_%.1f_%s.png' % (DLAT, sznnm), check_if_seeded=True)
      genhurr = plot_lat_binned(szn_dfs, bininfo, totyrs, genesis_pts_1d, 'genesis points of hurricanes', 'genhurr_pt_dens_%.1f_%s.png' % (DLAT, sznnm), peak_wspd_thresh=33)
      lys = plot_lat_binned(szn_dfs, bininfo, totyrs, lysis_pts_1d, 'lysis points', 'lys_pt_dens_%.1f_%s.png' % (DLAT, sznnm))
      us_lys = plot_lat_binned(szn_dfs, bininfo, totyrs, lysis_pts_1d, 'unseeded lysis points', 'us_lys_pt_dens_%.1f_%s.png' % (DLAT, sznnm), check_if_unseed=True)
      miss_us = plot_lat_binned(szn_dfs, bininfo, totyrs, missed_unseed, 'TCs missed by unseed', 'miss_us_pt_dens_%.1f_%s.png' % (DLAT, sznnm))

      #count storms by how many times they were unseeded
      us0 = plot_lat_binned(szn_dfs, bininfo, totyrs, lysis_pts_1d, 'stms with 0 unseeds', 'us0_%.1f_%s.png' % (DLAT, sznnm), unseed_attempts=0)
      us1 = plot_lat_binned(szn_dfs, bininfo, totyrs, lysis_pts_1d, 'stms with 1 unseeds', 'us1_%.1f_%s.png' % (DLAT, sznnm), unseed_attempts=1)
      us2 = plot_lat_binned(szn_dfs, bininfo, totyrs, lysis_pts_1d, 'stms with 2 unseeds', 'us2_%.1f_%s.png' % (DLAT, sznnm), unseed_attempts=2)
      us3 = plot_lat_binned(szn_dfs, bininfo, totyrs, lysis_pts_1d, 'stms with 3 unseeds', 'us3_%.1f_%s.png' % (DLAT, sznnm), unseed_attempts=3)
   
      outdss.append(xr.Dataset(data_vars=dict(uniq=uniq, h6all=h6all, h6hurr=h6hurr, h6maj=h6maj, ace=ace, gen=gen, sd_gen=sd_gen, genhurr=genhurr, lys=lys, us_lys=us_lys,\
                us0=us0, us1=us1, us2=us2, us3=us3), attrs=dict(dlat=DLAT)).expand_dims(season=[sznnm]))
   
      if PLOTSEEDS:
         evdss = []
         for ii, ev in enumerate(events):
            if ev is None:
               evdss.append(None)
               continue
            mydf, accum_ev = evdfs[ii], None
            if ev[1] == 'us':
               accum_ev = accum_bin_map_1d(mydf[mydf['dt'].dt.month.isin(sznmos)], bininfo, 'psmin', lambda ps: not np.isnan(ps), dtype=np.int_, rettype=np.int_)
            elif ev[1] == 'sd':
               accum_ev = accum_bin_map_1d(mydf[mydf['dt'].dt.month.isin(sznmos)], bininfo, 'lat', lambda lat: 1, dtype=np.int_, rettype=np.int_)

            evdens = xr.DataArray(accum_ev / (bininfo[2] / 1e6 / 1e3**2) / totyrs[ii], dims=['lat'], coords=[bininfo[1]])   
            outdss[-1] = outdss[-1].assign(variables={labels[ii] + ev[1]: evdens})
   
      xr.concat(outdss, dim='season').to_netcdf(os.path.join(DOUT, 'tcdens.nc'))

def plot_lat_binned(dfs, bininfo, totyrs, varfunc, ylabel, filo, norm='per million km2 per yr', **fkwargs):
   pltdat = [np.zeros_like(bininfo[1]) for _ in dfs]
   for ii, df in enumerate(dfs):
      for stm in df['stmnum'].unique():
         pltdat[ii] += varfunc(df[df['stmnum'] == stm], bininfo, **fkwargs)
      #print('\t', ylabel, '\t\t\t', pltdat[ii][nhmask].sum(), '\t\t\t', pltdat[ii][shmask].sum())
      with open(os.path.join(DOUT, 'sznl_climo.csv'), 'a') as pfil:
         print(df.index[0].month, ylabel, labels[ii], pltdat[ii][nhmask].sum(), pltdat[ii][shmask].sum(), totyrs[ii], sep=',', file=pfil)
      #print(pltdat[ii].shape)
      #print(pltdat[ii])
      #print(pltdat[ii].ravel()[nhmask.ravel()])
      pltdat[ii] = pltdat[ii] / (bininfo[2] / 1e6 / 1e3**2) / totyrs[ii]
      
      plt.plot(np.sin(np.deg2rad(bininfo[1])), pltdat[ii], label=labels[ii])
      plt.xticks(np.sin(np.deg2rad(LATLABS)), labels=LATLABS.astype(np.int_))
      plt.xlabel('lat')
      if varfunc == sixhrly_fixes_of_storm_1d:
         plt.ylim(*FIXYLIM)
      plt.ylabel('%s %s' % (ylabel, norm))
   plt.legend(fontsize=12)
   plt.savefig(os.path.join(DOUT, filo), bbox_inches='tight')
   plt.close()

   return xr.DataArray(np.array(pltdat), dims=['run', 'lat'], coords=dict(run=labels, lat=bininfo[1]), attrs=dict(norm=norm, long_name=ylabel))

#make an array of lat bin interfaces
#In: latitude bounds [deg], desired lat spacing [deg] at the equator
def make_bin_grid_1d(latbnds=(-50., 50.), dlat=0.5):
   sinbnds = (np.sin(np.deg2rad(latbnds[0])), np.sin(np.deg2rad(latbnds[1])))
   dsin = np.sin(np.deg2rad(dlat))
   true_dsin = (sinbnds[1] - sinbnds[0]) / ((sinbnds[1] - sinbnds[0]) // dsin + 1) #choose a spacing so that the interval is evenly divided

   bini = np.arange(sinbnds[0], sinbnds[1] + 1e-3, true_dsin) #bin interfaces in sin(lat) space
   binm = (bini[1:] + bini[:-1]) / 2
   area = 2 * np.pi * a_e**2 * true_dsin

   return np.rad2deg(np.arcsin(bini)), np.rad2deg(np.arcsin(binm)), area

#Binary search latitude bins and accumulate the desired quantity
#In: dataframe subset for 1 storm, bininfo (bin interfaces, bin midpoints, bin area [m^2]), column name, func to apply to column values, accum data type, return data type
def accum_bin_map_1d(stmdf, bininfo, col, afunc, dtype=np.float64, rettype=np.float64):
   outarr = np.zeros_like(bininfo[1], dtype=dtype)
   accvals = stmdf[col].apply(afunc)
   for ii in stmdf.index:
      if np.isnan(stmdf['lat'][ii]):
         continue
      iidx = np.searchsorted(bininfo[0], stmdf['lat'][ii])
      if iidx == 0 or iidx == bininfo[0].shape[0] - 1:
         continue
      outarr[iidx - 1] += accvals[ii]
   return outarr.astype(rettype)

#bool-to-int array of whether the storm hit the bin with midpoint
#i.e., return should be 0 or 1 divided by bin area
#In: dataframe subset for 1 storm, bininfo (bin interfaces, bin midpoints, bin area [m^2])
def unique_tracks_of_storm_1d(stmdf, bininfo):
   return accum_bin_map_1d(stmdf, bininfo, 'lat', lambda x: 1, dtype=np.bool_, rettype=np.int_)

def sixhrly_fixes_of_storm_1d(stmdf, bininfo, wspd_thresh=0):
   return accum_bin_map_1d(stmdf, bininfo, 'wspd', lambda wspd: int(wspd >= wspd_thresh), dtype=np.int_, rettype=np.int_)

def bin_ace_of_storm_1d(stmdf, bininfo):
   ace = lambda wspd: 1e-4 * (MS2KT * wspd)**2
   return accum_bin_map_1d(stmdf, bininfo, 'wspd', ace)

def genesis_pts_1d(stmdf, bininfo, peak_wspd_thresh=0, check_if_seeded=False):
   #gendf = stmdf.groupby('stmnum').first().reset_index() #wrong for seasonally split dataframes
   #return sixhrly_fixes_of_storm_1d(gendf, bininfo)
   mydf = stmdf
   if check_if_seeded:
      mydf = mydf[mydf['is_seeded'] == True] if 'is_seeded' in mydf.columns else mydf.iloc[0:0]
   return accum_bin_map_1d(mydf[mydf['max_lft_wspd'] >= peak_wspd_thresh], bininfo, 'isgen', lambda ig: int(ig), dtype=np.int_, rettype=np.int_)

def lysis_pts_1d(stmdf, bininfo, check_if_unseed=False, unseed_attempts=None):
   #lysdf = stmdf.groupby('stmnum').last().reset_index()
   #return sixhrly_fixes_of_storm_1d(lysdf, bininfo)
   mydf = stmdf
   if unseed_attempts is not None:
      if 'unseed_pt' in mydf.columns:
         #uscnt = mydf.groupby(mydf['stmnum'])['unseed_pt'].sum() #will only count unseed attempts strictly within season is stmdf is seasonal
         #print(uscnt, uscnt.index)
         #print(uscnt)
         #print(uscnt[uscnt == unseed_attempts].index)
         #mydf = mydf[mydf['stmnum'].isin(uscnt[uscnt == unseed_attempts].index)]
         #print(mydf['stmnum'])
         uscnt_per_row = mydf.groupby('stmnum')['unseed_pt'].transform('sum')
            
         # Filter rows where the parent storm matches the count
         mydf = mydf[uscnt_per_row == unseed_attempts]
      else:
         mydf = mydf.iloc[0:0]
   if check_if_unseed:
      mydf = mydf[mydf['unseed_pt'] == True] if 'unseed_pt' in mydf.columns else mydf.iloc[0:0]
   return accum_bin_map_1d(mydf, bininfo, 'islys', lambda il: int(il), dtype=np.int_, rettype=np.int_)

#count the number of storms that underwent lysis in this season and who were at some point eli
def missed_unseed(stmdf, bininfo):
   mydf = stmdf[stmdf['islys']]
   if 'us_elig' in mydf.columns:
      mydf = mydf[(mydf['istouched'] == 0) & (mydf['us_elig'])]
   else:
      mydf = mydf.iloc[0:0]
   return accum_bin_map_1d(mydf, bininfo, 'islys', lambda il: int(il), dtype=np.int_, rettype=np.int_)

if __name__ == '__main__':
   main()
