import xarray as xr
import uxarray as ux
import numpy as np
#from TEM_sznl import trapint
import matplotlib.pyplot as plt
import matplotlib.colors as colors

REST = '/glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1.cam.r.0005-02-19-00000.nc.ORIG.nc.plev'
CAMGRID = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

#compute the trapezoidal integral given integrand values and 1D coordinate values
#start at a zero integral value for the 1st coord value
def trapint(igrnd, coords):
   print(igrnd.name)
   dim = list(coords.coords)[0]
   #sz = coords.shape[0]
   igral = xr.DataArray(np.zeros(igrnd.shape, dtype=np.float64), dims=igrnd.dims, coords=igrnd.coords)
   for i, pt in enumerate(coords[1:]):
      prev = {dim: i} #use [] to index DataArrays
      curr = {dim: pt} #use .loc[] to index DataArrays
      step = pt - coords[prev]
      igral.loc[curr] = igral[prev] + (igrnd.loc[curr] + igrnd[prev]) / 2 * step
   return igral


ds = ux.open_mfdataset(CAMGRID, REST)
vzm = ds['V'].zonal_mean(lat=(-90, 90, .5))

coslat = np.cos(np.deg2rad(vzm.latitudes))
igrnd = vzm.isel(time=0) * 2 * np.pi * 6.371e6 * coslat / 9.81
igrnd = igrnd.assign_coords(coords=dict(plev=ds['plev']))
print(ds['plev'])
PSI_EM = trapint(igrnd, ds['plev'])

sinlat = np.sin(np.deg2rad(vzm.latitudes))
contourfkwargs = {'cmap': 'coolwarm', 'levels': 2.**np.arange(-2, 7, 1), 'norm': colors.SymLogNorm(2.**-1)}
contourfkwargs['levels'] = np.concatenate((-contourfkwargs['levels'][::-1], contourfkwargs['levels']))
contourkwargs = {'colors': 'black', 'levels': 2.**np.arange(-2, 7, 1)}
contourkwargs['levels'] = np.concatenate((-contourkwargs['levels'][::-1], contourkwargs['levels']))
csf = plt.contourf(sinlat, ds['plev'], PSI_EM / 1e10, **contourfkwargs)
plt.contour(sinlat, ds['plev'], PSI_EM / 1e10, **contourkwargs)
plt.yscale('log')
plt.gca().invert_yaxis()
plt.colorbar(csf)
plt.show()
