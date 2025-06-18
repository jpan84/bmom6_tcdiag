import numpy as np
import xarray as xr
import uxarray as ux

SZNS = ['DJF', 'MAM', 'JJA', 'SON']
#create months to seasons weight matrix
monlen = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.int_)
wmat = np.repeat(monlen[None, :], 4, axis=0)
for rr in range(wmat.shape[0]):
   for cc in range(wmat.shape[1]):
      mod = (cc + 1) % 12
      if mod >= 3 * (rr + 1) or mod < 3 * rr:
         wmat[rr, cc] = 0
wmat = wmat / wmat.sum(axis=1)[:, None]
#print(wmat)

def monthly2sznl(da, monnm='month'):
   da = da.transpose(monnm, ...)
   sznarr = wmat @ da.data

   coords = dict(da.coords)
   coords.pop(monnm)
   coords['season'] = SZNS

   outda = None
   if type(da) == xr.DataArray:
      outda = xr.DataArray(sznarr, coords=coords)
   if type(da) == ux.DataArray:
      outda = ux.UxDataArray(sznarr, coords=coords, uxgrid=da.uxgrid)

   return outda
