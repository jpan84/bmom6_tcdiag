#from h0a to area-averaged time series, allowing for lat selection and vertical integration
#keep everything in CESM output units

import os
import sys
from functools import reduce
import operator
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
from paths import ARCHRT, CAMGR
import consts as c
import numpy as np
import xarray as xr
import uxarray as ux

#DIRI = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250702_unseed2hPa6m/atm/hist/'
#DIRI = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm/hist/'
#DIRI = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist/'
#DIRI = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/atm/hist/'
#DIRI = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist/'

HTAPE = '*.h0a.*.nc'

#LATBNDS = (-90, 90)
#LATBNDS = (5, 35)
LATBNDS = (10, 30)
LATBNDS = (-30, -10)

RAWV2D = ['PS', 'TGCLDIWP', 'TGCLDLWP', 'QFLX', 'TS', 'FLNT', 'FSNT', 'LWCF', 'SWCF', 'FLNS', 'FSNS']
UDEF2D = dict(TAUAM=[(c.a_e, 'coslat', 'TAUX')], SHU=[(1, 'FSNS'), (-1, 'FLNS'), (-1, 'SHFLX'), (-1, 'LHFLX')],
              AHSRC=[(1, 'FSNT'), (-1, 'FSNS'), (1, 'FLNS'), (-1, 'FLNT'), (1, 'SHFLX'), (c.lv, c.rho_w, 'PRECC'), (c.lv, c.rho_w, 'PRECL')],
              LE = [(c.lv, 'TMQ')], PRECT=[(1, 'PRECC'), (1, 'PRECL')])
RAWV3D = None
UDEF3D = dict(DSE=[(c.cp, 'T'), (c.g, 'Z3')], LE_3D=[(c.lv, 'Q')], KE=[(0.5, 'UU'), (0.5, 'VV')], KE_MEAN=[(0.5, 'U', 'U'), (0.5, 'V', 'V')],
              AM=[(c.a_e, 'coslat', 'U')], CLD_FT=[(1, 'CLDICE'), (1, 'CLDLIQ')], CLD_BL=[(1, 'CLDICE'), (1, 'CLDLIQ')])
PTROP = (1e4, 1.1e5)
PBNDS = dict(DSE=PTROP, LE_3D=PTROP, KE=PTROP, KE_MEAN=PTROP, CLD_FT=(1e4, 7e4), CLD_BL=(7e4, 1.1e5))

def main(diri):
   print('Opening dataset', diri, '...')
   ds = ux.open_mfdataset(CAMGR, os.path.join(diri, HTAPE))

   print('Setting up coords...')
   aterm = ds['hyai'] * c.P0
   bterm = ds['hybi'] * ds['PS']
   p_ilev = aterm + bterm
   dp3d = p_ilev.diff('ilev').rename(dict(ilev='lev')).assign_coords(lev=ds['lev'])
   ds = ds.assign(variables=dict(dp3d=dp3d, p_ilev=p_ilev, coslat=np.cos(np.deg2rad(ds['lat']))))

   print('Deriving 2D quantities...')
   qtys_2d = {vn: ds[vn] for vn in RAWV2D}
   qtys_2d = qtys_2d | {kk: apply_udef_template(ds, vv) for kk, vv in UDEF2D.items()}

   print('Vertically integrating 3D fields...')
   qtys_3d = {kk: apply_udef_template(ds, vv) for kk, vv in UDEF3D.items()}
   dp_integ = {}
   for kk, vv in qtys_3d.items():
      mydp = ds['dp3d']
      if kk in PBNDS:
         mydp = get_dp_clipped(ds, *PBNDS[kk])
      dp_integ[kk] = (vv * mydp).sum(dim='lev') / c.g

   print('Area averaging...')
   aavg = {}
   for kk, vv in (qtys_2d | dp_integ).items():
      latsel = vv.subset.constant_latitude_interval(LATBNDS)
      aavg[kk] = latsel.weighted_mean()
      print('\t', type(aavg[kk]))

   print('Saving output...')
   outds = xr.Dataset(data_vars=aavg)
   outds.to_netcdf(os.path.join(diri, '../AAVG_%.1f_%.1f_rev1.nc' % (LATBNDS)))

   print(__file__, 'done.')

def apply_udef_template(ds, template):
   def get_product(ds, tup):
      vals = [ds[v] if isinstance(v, str) else v for v in tup]
      return reduce(operator.mul, vals)

   #Calculate the raw product (your current logic)
   terms = [get_product(ds, t) for t in template]
   total = sum(terms)
   
   return total

def get_dp_clipped(ds, pt, pb):
   pi_clip = ds['p_ilev'].where(ds['p_ilev'] < pb, pb).where(ds['p_ilev'] > pt, pt)
   return pi_clip.diff('ilev').rename(dict(ilev='lev')).assign_coords(lev=ds['lev'])

if __name__ == '__main__':
   for ar in ARCHRT:
      main(os.path.join(ar, 'atm/hist'))
