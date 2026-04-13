#from h0a to area-averaged time series, allowing for lat selection and vertical integration
#keep everything in CESM output units

import operator
import consts as c
import numpy as np
import xarray as xr
import uxarray as ux

LATBNDS = (5, 35)

RAWV2D = ['PS', 'TGCLDIWP', 'TGCLDLWP', 'PRECT', 'QFLX', 'TS', 'FLNT', 'FSNT', 'LWCF', 'SWCF', 'FLNS', 'FSNS']
UDEF2D = dict(TAUAM=[(c.a, 'coslat', 'TAUX')], SHU=[(1, 'FSNS'), (-1, 'FLNS'), (-1, 'SHFLX'), (-1, 'LHFLX')],
              AHSRC=[(1, 'FSNT'), (-1, 'FSNS'), (1, 'FLNS'), (-1, 'FLNT'), (1, 'SHFLX'), (c.lv, c.rho_w, 'PRECT')],
              LE = [(c.lv, 'TMQ')])
RAWV3D = None
UDEF3D = dict(DSE=[(c.cp, 'T'), (c.g, 'Z3')], KE=[(1, 'UU'), (1, 'VV')], KE_MEAN=[(1, 'U', 'U'), (1, 'V', 'V')],
              AM=[(c.a, 'coslat', 'U')], CLD_FT=[(1, 'CLDICE'), (1, 'CLDLIQ')], CLD_BL=[(1, 'CLDICE'), (1, 'CLDLIQ')])
PBNDS = dict(CLD_FT=(1e4, 7e4), CLD_BL=(7e4, 1.1e5))

def main():
   print('Opening dataset...')
   ds = ux.open_mfdataset(os.path.join(DIRI, HTAPE))

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
   for kk, vv in qtys_2d | dp_integ:
      latsel = vv.subset.constant_latitude_interval(LATBNDS)
      aavg[kk] = latsel.weighted_mean()

   print('Saving output...')
   outds = xr.Dataset(data_vars=aavg)
   outds.to_netcdf(os.path.join(DIRI, 'AAVG_%.1f_%.1f.nc' % (LATBNDS)))

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
   main()

#class var2d:
#   def __init__(self, list: template, varnm=None, pbnds=None):
#class h0a_reduce:
#   def __init__(self, ds):
#      self.ds = ds
#      #TODO: add dp3d and coslat
#      aterm = self.ds['hyai'] * c.P0
#      bterm = self.ds['hybi'] * self.ds['PS']
#      p_ilev = aterm + bterm
#      dp3d = p_ilev.diff('ilev').rename(dict(ilev='lev'))
#
#      self.ds = self.ds.assign(variables=dict(dp3d=dp3d, p_ilev=p_ilev, coslat=np.cos(np.deg2rad(self.ds['lat']))))
#
#   def apply_udef_template(self, template=None, varnm=None):
#      # 1. Calculate the raw product (your current logic)
#      terms = [self._get_product(t) for t in definition[:-1]]
#      total = sum(terms)
#      
#      # It's a 2D variable (like PS or SST): return as is
#      return total
#
#   def _get_product(self, tup):
#      vals = [self.ds[v] if isinstance(v, str) else v for v in tup]
#      return reduce(operator.mul, vals)
#
#   def dp_integ(self, pbnds=None):
#      # 2. Automated Vertical Handling
#      if 'lev' in total.dims or 'ilev' in total.dims:
#          # It's a 3D variable: it needs vertical integration (dp/g)
#          # This is where your conservative dp logic lives!
#          return (total * self.ds['dp'] / self.g).sum(dim='lev')
#      
