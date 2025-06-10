#Joshua Pan orig from Oct 2023 for docn runs
#Updated May 2025 to use uxarray non-conservative zonal mean
#plot and compare zonal-mean time-averaged quantities

#TOPDIR = '/glade/scratch/zarzycki/archive/'
#DIR1 = 'QPC5-ne30np4-aquap10-seed3x3/atm/hist'
#DIR2 = 'QPC5-ne30np4-aquap10-unseed/atm/hist'
#FN = 'QPC5-ne30np4-aquap10-*.cam.h0.*regrid.nc'
OUTDIR = 'linevslat_250415_unseed_minus_ctrl/'
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
import pltsettings

### hist file params
MODE = 'CAM'
ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed'
HISTS = 'atm/hist/'
H0 = '*.cam.h0a.*.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

'''
MODE = 'MOM'
HISTS = 'ocn/hist/'
H0 = '*mom6.hm*[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'
'''

DO_DIFF = True
CASE2 = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl'

grpby = 'season' #month
SZNS = ['DJF', 'MAM', 'JJA', 'SON']
zmlats = (-90, 90, 0.5)
LATLAB = np.array([-90., -60., -30., 0., 30., 60., 90.])
lncolors = plt.cm.jet(np.linspace(0, 1, 12 if grpby == 'month' else 4 if grpby == 'season' else None))
#TODO: allow diffing between cases and selecting of months/seasons
SKIP = {'AEROD_v'}
USER_DEF = {'RESTOM', 'PRECT', 'NCF'}

HEMISYM = {'FSNS', 'FLNS', 'LHFLX', 'SHFLX'}
DO_SYM = False #only works for season, not month

#create months to seasons weight matrix
mlnl = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.int_)
wmat = np.repeat(mlnl[None, :], 4, axis=0)
for rr in range(wmat.shape[0]):
   for cc in range(wmat.shape[1]):
      mod = (cc + 1) % 12
      if mod >= 3 * (rr + 1) or mod < 3 * rr:
         wmat[rr, cc] = 0
wmat = wmat / wmat.sum(axis=1)[:, None]
print(wmat)

def main():
   if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

   pltsettings.set1()

   print('Opening datasets...')
   
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
      if str(dv) in SKIP or DO_SYM and str(dv) not in HEMISYM:
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

if __name__ == '__main__':
   main()
