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
OUTDIR = 'linevslat_new_sznl_diff'
ARCHV = '/glade/derecho/scratch/jpan/archive/'
HISTS = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/*.h0a.*.nc'
CASES = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']
ALIASES = ['UNSEED', 'CTRL', 'SEED']
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

DO_DIFF = True

zmlats = (-90, 90, 0.5)
LATLAB = np.array([-90., -60., -30., 0., 30., 60., 90.])
lncolors = ['blue', 'orange']
#TODO: allow diffing between cases and selecting of months/seasons
SKIP = {'AEROD_v'}
USER_DEF = {'RESTOM', 'PRECT', 'NCF'}

#h1i mode
#HISTS = '/glade/derecho/scratch/jpan/archive/%s/atm/hist/*.h1i.0010*.nc'
H1I = False
H1IVARS = {'CF500'}
#USER_DEF = {'CF500'}

def main():
   if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

   pltsettings.set1()

   print('Opening datasets...')
   dss = [ux.open_mfdataset(camgrid, HISTS % cs) for cs in CASES]
   dvset = set([str(dv) for dv in dss[0].data_vars])
   dvset = dvset | USER_DEF
   outds = None
   for dv in dvset:
      if str(dv) in SKIP or H1I and str(dv) not in H1IVARS:
         print('Skipping %s...' % dv)
         continue
      
      #get/compute vars of shape (time, ncol)
      das = None
      if str(dv) in USER_DEF:
         das = [udef(ds, dv) for ds in dss]
      elif set(dss[0][dv].dims) == HISTDIMS:
         das = [ds[dv] for ds in dss]
      else:
         print('Skipping %s...' % dv)
         continue
      print('Processing %s...' % dv)

      print('\tComputing means...')
      monmeans = [da.groupby('time.month').mean() for da in das]
      monzm = [da.zonal_mean(lat=zmlats) for da in monmeans] #monthly zonal means
      #print(monzm[0])
      #monzm = [da.assign_coords(month=monmeans[0]['month'].data) for da in monmeans] #fix uxarray zonal_mean() dropping coords
      #print(monzm[0])
      sznzm = [sznl_funcs.monthly2sznl(da) for da in monzm] #shape (season, ncol)
      sznzm = [sznl_funcs.stack_hemi_sznl(da) for da in sznzm]
      sinlat = np.sin(np.deg2rad(sznzm[0]['latitudes']))

      if DO_DIFF:
         sznzm[0] = sznzm[0] - sznzm[1]
         sznzm[2] = sznzm[2] - sznzm[1]

      print('\tPlotting...')
      plt.rcParams['figure.figsize'] = (30, 6)
      subplot_kw = dict(xlim=(-1, 1), sharey=(not DO_DIFF))
      fig, axes = plt.subplots(1, 3, subplot_kw=subplot_kw)

      for ii, pltda in enumerate(sznzm):
         ax = axes[ii]
         for tt, szn in enumerate(pltda['season']):
            ax.plot(sinlat, pltda.sel(season=szn), label=str(szn.values), color=lncolors[tt])
         #plt.plot(sinlat, line1.values)
         #if str(dv) == 'TS' and not DO_DIFF:
         #   plt.hlines(273.15 + 26.5, -1, 1, colors='red', linestyles='dashed')
         if DO_DIFF and (ii == 0 or ii == 2):
            ax.hlines(0, -1, 1, colors='black', linestyles='dashed')
         ax.set_xlabel('Lat [°]')
         try:
            ax.set_ylabel(str(dv) + ' [' + str(das[1].units) + ']')
         except:
            ax.set_ylabel(dv)
         ax.set_title(ALIASES[ii] + (' difference' if (DO_DIFF and (ii == 0 or ii == 2)) else '') )
         ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=12))
         ax.set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))

      if DO_DIFF:
         ###make the difference panels share y
         #ylims = [ax.get_ylim() for ax in axes]
         #ylims.pop(1)
         #ylims = np.array(ylims)
         #miny, maxy = ylims.min(), ylims.max()
         #axes[0].set_ylim(miny, maxy)
         #axes[2].set_ylim(miny, maxy)
         for aa in [axes[0], axes[2]]:
            ylims = aa.get_ylim()
            maxy = max(np.abs(ylims))
            aa.set_ylim(-maxy, maxy)
      plt.savefig(os.path.join(OUTDIR, '%s_sznl.png' % dv), bbox_inches='tight')
      plt.close()

      #Add to output dataset/netcdf
      if outds is None:
         outds = xr.Dataset(data_vars={dv: xr.concat(sznzm, dim=xr.Variable('case', ALIASES))})
         outds.attrs['difference'] = str(DO_DIFF)
      else:
         outds = outds.assign(variables={dv: xr.concat(sznzm, dim=xr.Variable('case', ALIASES))})

   print('Saving to .nc...')
   outds.to_netcdf(os.path.join(OUTDIR, 'sznlzm.nc'))

   print(sys.argv[0], 'done.')


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
