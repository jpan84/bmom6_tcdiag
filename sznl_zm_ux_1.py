#Joshua Pan orig from Oct 2023 for docn runs
#Updated May 2025 to use uxarray non-conservative zonal mean
#plot and compare zonal-mean time-averaged quantities

HISTDIMS = set(['time', 'n_face']) #cam-SE

import os
import sys
import cftime
import datetime as dt
import numpy as np
import uxarray as ux
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import sznl_funcs
import pltsettings

### hist file params
OUTDIR = 'linevslat_new_sznl'
ARCHV = '/glade/derecho/scratch/jpan/archive/'
HISTS = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/*.h0a.*.nc'
CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
ALIASES = ['UNSEED', 'CTRL', 'SEED']
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

DO_DIFF = False

zmlats = (-90, 90, 0.5)
LATLAB = np.array([-90., -60., -30., 0., 30., 60., 90.])
lncolors = plt.cm.jet(np.linspace(0, 1, 4))
#TODO: allow diffing between cases and selecting of months/seasons
SKIP = {'AEROD_v'}
USER_DEF = {'RESTOM', 'PRECT', 'NCF'}

def main():
   if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

   pltsettings.set1()

   print('Opening datasets...')
   dss = [ux.open_mfdataset(camgrid, HISTS % cs) for cs in CASES]
   dvset = set([str(dv) for dv in dss[0].data_vars])
   dvset = dvset | USER_DEF
   for dv in dvset:
      if str(dv) in SKIP:
         print('Skipping %s...' % dv)
         continue

      print('Processing %s...' % dv)
      #get/compute vars of shape (time, ncol)
      das = None
      if str(dv) in USER_DEF:
         das = [udef(ds, dv) for ds in dss]
      elif set(dss[0][dv].dims) == HISTDIMS:
         das = [ds[dv] for ds in dss]

      print('     Computing means...')
      monmeans = [da.groupby('time.month').mean() for da in das]
      monzm = [tmeans.zonal_mean(lat=zmlats).to_xarray()] #monthly zonal means
      sznzm = [sznl_funcs.monthly2sznl(da) for da in monzm]
      sinlat = np.sin(np.deg2rad(sznzm[0]['latitudes']))

      print('     Plotting...')
      plt.rcParams['figure.figsize'] = (30, 6)
      subplot_kw = dict(xlim=(-1, 1), sharey=(not DO_DIFF))
      fig, axes = plt.subplots(1, 3, subplot_kw=subplot_kw)

      for ii, pltda in enumerate(sznzm):
         ax = axes[ii]
         for tt, szn in enumerate(pltda['season']):
            ax.plot(sinlat, pltda.sel(season=szn), label=szn, color=lncolors[tt])
         #plt.plot(sinlat, line1.values)
         #if str(dv) == 'TS' and not DO_DIFF:
         #   plt.hlines(273.15 + 26.5, -1, 1, colors='red', linestyles='dashed')
         #if DO_DIFF:
         #   plt.hlines(0, -1, 1, colors='black', linestyles='dashed')
         ax.set_xlabel('Lat [°]')
         ax.set_ylabel(dv)
         ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=12))
         ax.set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
      plt.savefig(os.path.join(OUTDIR, '%s_sznl.png' % dv), bbox_inches='tight')
      plt.close()

   print(sys.argv[0], 'done.')

   exit()
   ##################################################################
   
   ds1 = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HISTS, H0))
   ds2 = None
   if DO_DIFF:
      ds2 = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE2, HISTS, H0))
   #ds1 = ds1.assign_coords(coords=dict(time=ds1.time - dt.timedelta(days=1))) #timestamp of e.g., 0001-01.nc is Feb 1, so subtract a month from time coord
   #ds1 = ds1.assign(dict(U200=ds1.U.sel(lev=200, method='nearest')))
   #!OCEAN ds1 = ds1.isel(z_l=0)
   print(ds1)

   dvset = set([str(dv) for dv in ds1.data_vars])
   dvset = dvset | USER_DEF
   for dv in dvset:
      if str(dv) in SKIP or DO_SYM and str(dv) not in HEMISYM or str(dv) not in {'PRECT', 'OMEGA500', 'TC8GCD5MPS', 'CF500'}: #!HIFREQ
         print('Skipping %s...' % dv)
         continue
      if str(dv) in USER_DEF or set(ds1[dv].dims) == HISTDIMS:
         print('Plotting %s...' % dv)
         da = None
         if dv in USER_DEF:
            da = udef(ds1, dv)
         else:
            da = ds1[dv]
         if DO_DIFF:
            if dv in USER_DEF:
               da = da - udef(ds2, dv)
            else:
               da = da - ds2[dv]
         tmeans = da.groupby('time.month').mean()
         #lines = None
         #if
         lines = tmeans.zonal_mean(lat=zmlats)
         sinlat = np.sin(np.deg2rad(lines.latitudes))
         if grpby == 'season':
            print(lines.shape)
            lines = wmat @ lines.values
         for tt in range(lines.shape[0]):
            plt.plot(sinlat, lines[tt], label=SZNS[tt] if grpby == 'season' else tt, color=lncolors[tt])
         #plt.plot(sinlat, line1.values)
         if str(dv) == 'TS' and not DO_DIFF:
            plt.hlines(273.15 + 26.5, -1, 1, colors='red', linestyles='dashed')
         if DO_DIFF:
            plt.hlines(0, -1, 1, colors='black', linestyles='dashed')
         if str(dv) == 'PRECT': #!HIFREQ
            plt.ylim(0, 2.2e-7)
         if str(dv) == 'TC8GCD5MPS':
            plt.ylim(0, 0.14)
         plt.xlabel('Lat [°]')
         plt.ylabel(dv)
         plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=12))
         plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
         plt.savefig(os.path.join(OUTDIR, '%s_sznl.png' % dv), bbox_inches='tight')
         #plt.tight_layout()
         #plt.show()
         plt.close()

         if DO_SYM:
            nhstart = int(lines.shape[1] // 2) + 1 #odd number of lats, sym about Eq
            shend = nhstart - 1
            nh = lines[:, nhstart:]
            sh = np.roll(lines[:, :shend], 2, axis=0)[:, ::-1] #align seasons and reflect across Eq
            hemidiff = nh - sh
            for tt in range(hemidiff.shape[0]):
               plt.plot(sinlat[nhstart:], hemidiff[tt], label='NH %s' % SZNS[tt], color=lncolors[tt])
            plt.hlines(0, -1, 1, colors='black', linestyles='dashed')
            plt.vlines(np.sin(np.deg2rad(73.75)), 0, 10, colors='gray', linestyles='dashed')
            plt.xlabel('Lat [°]')
            plt.ylabel(dv)
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=12))
            plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
            plt.xlim(0, 1)
            plt.savefig(os.path.join(OUTDIR, '%s_hemidiff.png' % dv), bbox_inches='tight')
            plt.close()
      else:
         print('Skipping %s...' % dv)

   print('%s done.' % sys.argv[0])

def udef(ds, dv):
   if dv == 'PRECT':
      return ds['PRECC'] + ds['PRECL']
   if dv == 'RESTOM':
      return ds['FSNT'] - ds['FLNT']
   if dv == 'NCF':
      return ds['LWCF'] + ds['SWCF']
   if dv == 'CF500':
      return ux.UxDataArray((ds['OMEGA500'] < 0).astype(np.float64).data, dims=ds['OMEGA500'].dims,\
             coords=ds['OMEGA500'].coords, uxgrid=ds.uxgrid)

if __name__ == '__main__':
   main()
