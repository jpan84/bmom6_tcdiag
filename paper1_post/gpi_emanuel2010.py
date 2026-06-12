import uxarray as ux
from uxarray.core import gradient
import xarray as xr
import numpy as np

from tcpyPI import pi as tcpi
from entropy_deficit_wcy import entropy_deficit
from metpy.calc import relative_humidity_from_specific_humidity as q2rh
from metpy.units import units

from paths import CAMGR
from consts import TCK, OM
import matplotlib.pyplot as plt

TSTFIL = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm/hist_onpres_gnupar/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.cam.h0a.0011-08.onpres.nc'

pintp = lambda x, p: x.interp(plev=p).drop('plev')

def main():
   ds = ux.open_mfdataset(CAMGR, TSTFIL)
   ds = ds.sortby('plev', ascending=False)
   ds = ds.chunk({'time': 1, 'plev': -1, 'n_face': -1})
   print(ds)

   SSTC = ds['TS'] - TCK
   MSL = ds['PSL'] / 100
   phPa = ds['plev'] / 100
   TC = ds['T'] - TCK
   rwv = ds['Q'] / (1 - ds['Q']) * 1000

   #PI vmax
   print('PI task graph...')
   vmax, _, ifl, _, _ = pi_xr_ux(SSTC, MSL, phPa, TC, rwv)
   vmax = ux.UxDataArray(vmax, uxgrid=ds.uxgrid)
   vmax_zm = vmax.zonal_mean(lat=(-50, 50, 2)).squeeze()
   #print('IFL', ifl.values)
   #print(vmax_zm.values)

   #plt.plot(vmax_zm.latitudes, vmax_zm)
   #plt.show()

   #Emanuel entropy deficit
   print('Chi task graph...')
   rhb = q2rh(ds['PS'], ds['TREFHT'], ds['QREFHT'])
   pm = 6e4 #Pa, mid-tropo reference level
   Tm = pintp(ds['T'], pm)
   rhm = q2rh(pm * units('Pa'), Tm, pintp(ds['Q'], pm))
   #print(rhb.max().values)
   #print(rhm.max().values)
   chi_ent = entropy_deficit(ds['TS'], ds['PSL'], ds['TREFHT'], rhb, pm,\
              Tm, rhm)
   #print(chi_ent.values)

   #850 absolute vort
   print('Vort task graph...')
   #rv = pintp(ds['U'], 8.5e4).curl(pintp(ds['V'], 8.5e4))
   rv = apply_ux_curl(pintp(ds['U'], 8.5e4), pintp(ds['V'], 8.5e4))
   #rv = ux.UxDataArray(rv, uxgrid=ds.uxgrid)
   #print(rv.max().values)
   va = (rv + 2 * OM * np.sin(np.deg2rad(ds['lat']))).clip(-5e-5, 5e-5)

   #250-850 shear magnitude
   print('Shear task graph...')
   ushr = pintp(ds['U'], 2e4) - pintp(ds['U'], 8.5e4)
   vshr = pintp(ds['V'], 2e4) - pintp(ds['V'], 8.5e4)
   shrmag = np.sqrt(ushr ** 2 + vshr ** 2)

   print('GPI task graph...')
   f1 = np.abs(va) ** 3
   f2 = chi_ent ** (-4 / 3)
   f3 = np.maximum(vmax, 0) ** 2
   f4 = (25 + shrmag) ** -4
   gpi = f1 * f2 * f3 * f4
   gpi = ux.UxDataArray(gpi, uxgrid=ds.uxgrid)
   gpi_zm = gpi.zonal_mean(lat=(-50, 50, 2)).squeeze()

   print('Computing and plotting...')
   plt.plot(gpi_zm.latitudes, gpi_zm)
   plt.show()


def pi_xr_ux(*piargs):
   return xr.apply_ufunc(tcpi, *piargs,\
            input_core_dims=[[], [], ['plev'], ['plev'], ['plev']],\
            output_core_dims=[[] for _ in range(len(piargs))],\
            vectorize=True, dask='parallelized',\
            output_dtypes=[float, float, int, float, float])

#TODO: try rewrapping np arrays uu and vv into UxDataArray
def ux_curl_oop(uu, vv):
   return uu.curl(vv)

#uu, vv should have the same coords (time, n_face)
#arrays must not be chunked along dim n_face
def apply_ux_curl(uu, vv, uxgrid):
   face_coords = np.array(
       [uxgrid.face_x.values, uxgrid.face_y.values, uxgrid.face_z.values]
   ).T
   face_lat = uxgrid.face_lat.values
   face_lon = uxgrid.face_lon.values

   node_coords = np.array(
       [uxgrid.node_x.values, uxgrid.node_y.values, uxgrid.node_z.values]
   )

   face_lon_rad = np.deg2rad(face_lon)
   face_lat_rad = np.deg2rad(face_lat)
   normal_lat = np.array(
       [
           -np.cos(face_lon_rad) * np.sin(face_lat_rad),
           -np.sin(face_lon_rad) * np.sin(face_lat_rad),
           np.cos(face_lat_rad),
       ]
   ).T
   normal_lon = np.array(
       [
           -np.sin(face_lon_rad),
           np.cos(face_lon_rad),
           np.zeros_like(face_lon_rad),
       ]
   ).T

   coordargs = (uxgrid.n_face, face_coords, uxgrid.edge_face_connectivity.values,\
                uxgrid.face_node_connectivity.values, uxgrid.node_edge_connectivity.values,\
                face_lat, face_lon, node_coords, normal_lon, normal_lat)

   callgrad = lambda da: xr.apply_ufunc(gradient._compute_gradients_on_faces, da, *coordargs,\
                 input_core_dims=[['n_face']] + [[] for _ in range(len(coordargs))],\
                 output_core_dims=[['n_face'], ['n_face']], vectorize=True,\
                 dask='parallelized', output_dtypes=[da.dtype, da.dtype])

   _, u_y = callgrad(uu)
   v_x, _ = callgrad(vv)

   return ux.UxDataArray(v_x - u_y, uxgrid=uxgrid)


if __name__ == '__main__':
   main()
