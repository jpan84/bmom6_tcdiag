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

   dims = list(da.dims)
   coords = dict(da.coords)
   print(da)
   print(dims)
   print(coords)
   if monnm in coords:
      coords.pop(monnm)
      dims.remove(monnm)
   dims = ['season'] + dims
   coords['season'] = SZNS

   outda = None
   if type(da) == xr.DataArray:
      outda = xr.DataArray(sznarr, dims=dims, coords=coords)
   if type(da) == ux.UxDataArray:
      outda = ux.UxDataArray(sznarr, dims=dims, coords=coords, uxgrid=da.uxgrid)

   return outda

#reduce a DataArray from 4 seasons to 2
#by mirroring across the equator and averaging
#lats must be increasing and symmetric about Eq (0 deg optional)
def stack_hemi_sznl(sznlda, antisym=False, sznnm='season', latnm='latitudes'):
   nhwarm = sznlda.sel({sznnm: ['JJA', 'SON']}) 
   shwarm = sznlda.sel({sznnm: ['DJF', 'MAM']})

   latidx = sznlda.dims.index(latnm)
   latflip = np.flip(shwarm.data, axis=latidx) * (-1 if antisym else 1)

   res = (nhwarm + latflip) / 2

   ###this code block is more like hemi difference logic
   #lats = sznlda[latnm]
   #nhstart = lats.size // 2
   #shend = nhstart + (lats.size % 2)

   #nhda = sznlda.isel({latnm: slice(nhstart, None)})
   #shda = sznlda.isel({latnm: slice(0, shend)})

   #sznidx = sznlda.dims.index(sznnm)
   #latidx = sznlda.dims.index(latnm)

   #sznshft = np.roll(shda.data, 2, axis=sznidx)
   #latflip = np.flip(sznshft, axis=latidx)

   return res
