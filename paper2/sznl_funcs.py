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
   #sznarr = wmat @ da.data
   sznarr = np.einsum('ij,j...->i...', wmat, da.data)

   dims = list(da.dims)
   coords = dict(da.coords)
   #print(da)
   #print(dims)
   #print(coords)
   if monnm in coords:
      coords.pop(monnm)
   if monnm in dims:
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
#lats must be symmetric about Eq (0 deg optional)
def stack_hemi_sznl(sznlda, antisym=False, sznnm='season', latnm='latitudes'):
   nhwarm = sznlda.sel({sznnm: ['JJA', 'SON']}) 
   shwarm = sznlda.sel({sznnm: ['DJF', 'MAM']})

   latidx = sznlda.dims.index(latnm)
   latflip = np.flip(shwarm.data, axis=latidx) * (-1 if antisym else 1)

   res = (nhwarm + latflip) / 2

   return res

def agg_time(da, antisym=False, sznnm='season', latnm='latitudes'):
   ymonmean = da.groupby('time.month').mean()
   twoszns = stack_hemi_sznl(monthly2sznl(ymonmean), antisym=antisym, latnm=latnm)
   halfyr = twoszns.mean(dim=sznnm)
   return halfyr
