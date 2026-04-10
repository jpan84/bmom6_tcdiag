#from h0a to area-averaged time series, allowing for lat selection and vertical integration
#keep everything in CESM output units

import consts as c
import uxarray as ux

#class var2d:
#   def __init__(self, list: template, varnm=None, pbnds=None):
class h0a_reduce:
   def __init__(self, ds):
      self.ds = ds

   def compute(self, list: template, pbnds=None):
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
