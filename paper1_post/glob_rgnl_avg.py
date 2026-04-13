#from h0a to area-averaged time series, allowing for lat selection and vertical integration
#keep everything in CESM output units

import operator
import consts as c
import uxarray as ux

RAWV2D = ['PS', 'TGCLDIWP', 'TGCLDLWP', 'PRECT', 'QFLX', 'TS', 'FLNT', 'FSNT', 'LWCF', 'SWCF', 'FLNS', 'FSNS']
UDEF2D = dict(TAUAM=[(c.a, 'coslat', 'TAUX')], SHU=[(1, 'FSNS'), (-1, 'FLNS'), (-1, 'SHFLX'), (-1, 'LHFLX')],
              AHSRC=[(1, 'FSNT'), (-1, 'FSNS'), (1, 'FLNS'), (-1, 'FLNT'), (1, 'SHFLX'), (c.lv, c.rho_w, 'PRECT')],
              LE = [(c.lv, 'TMQ')])
RAWV3D = None
UDEF3D = dict(DSE=[(c.cp, 'T'), (c.g, 'Z3')], KE=[(1, 'UU'), (1, 'VV')], KE_MEAN=[(1, 'U', 'U'), (1, 'V', 'V')],
              AM=[(c.a, 'coslat', 'U')], CLD_FT=[(1, 'CLDICE'), (1, 'CLDLIQ')], CLD_BL=[(1, 'CLDICE'), (1, 'CLDLIQ')])
PBNDS = dict(CLD_FT=(1e4, 7e4), CLD_BL=(7e4, 1.1e5))

#class var2d:
#   def __init__(self, list: template, varnm=None, pbnds=None):
class h0a_reduce:
   def __init__(self, ds):
      self.ds = ds
      #TODO: add dp3d and coslat

   def compute(self, template=None, varnm=None, pbnds=None):
      # 1. Calculate the raw product (your current logic)
      terms = [self._get_product(t) for t in definition[:-1]]
      total = sum(terms)
      
      # 2. Automated Vertical Handling
      if 'lev' in total.dims or 'ilev' in total.dims:
          # It's a 3D variable: it needs vertical integration (dp/g)
          # This is where your conservative dp logic lives!
          return (total * self.ds['dp'] / self.g).sum(dim='lev')
      
      # It's a 2D variable (like PS or SST): return as is
      return total

   def _get_product(self, tup):
      vals = [self.ds[v] if isinstance(v, str) else v for v in tup]
      return reduce(operator.mul, vals)
