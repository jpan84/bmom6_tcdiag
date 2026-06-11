import uxarray as ux
import xarray as xr

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
   ds = ds.chunk({'plev': -1, 'n_face': 777600 / 6 / 4})
   print(ds)

   SSTC = ds['TS'] - TCK
   MSL = ds['PSL'] / 100
   phPa = ds['plev'] / 100
   TC = ds['T'] - TCK
   rwv = ds['Q'] / (1 - ds['Q']) * 1000

   #PI vmax
   vmax, _, ifl, _, _ = pi_xr_ux(SSTC, MSL, phPa, TC, rwv)
   vmax = ux.UxDataArray(vmax, uxgrid=ds.uxgrid)
   vmax_zm = vmax.zonal_mean(lat=(-50, 50, 2)).squeeze()
   print('IFL', ifl.values)
   print(vmax_zm.values)

   #plt.plot(vmax_zm.latitudes, vmax_zm)
   #plt.show()

   #Emanuel entropy deficit
   rhb = q2rh(ds['PS'], ds['TREFHT'], ds['QREFHT'])
   pm = 6e4 #Pa, mid-tropo reference level
   Tm = pintp(ds['T'], pm)
   rhm = q2rh(pm * units('Pa'), Tm, pintp(ds['Q'], pm))
   print(rhb.max().values)
   print(rhm.max().values)
   chi_ent = entropy_deficit(ds['TS'], ds['PSL'], ds['TREFHT'], rhb, pm,\
              Tm, rhm)
   print(chi_ent.values)

   #850 absolute vort
   rv = pintp(ds['U'], 8.5e4).curl(pintp(ds['V'], 8.5e4))
   va = rv + 2 * OM * np.sin(np.deg2rad(ds['lat'])).clip(-5e-5, 5e-5)

   #250-850 shear magnitude
   ushr = pintp(ds['U'], 2e4) - pintp(ds['U'], 8.5e4)
   vshr = pintp(ds['V'], 2e4) - pintp(ds['V'], 8.5e4)
   shrmag = np.sqrt(ushr ** 2 + vshr ** 2)

   f1 = np.abs(va) ** 3
   f2 = chi_ent ** (-4 / 3)
   f3 = np.maximum(vmax, 0) ** 2
   f4 = (25 + shrmag) ** -4
   gpi = f1 * f2 * f3 * f4
   gpi_zm = gpi.zonal_mean(lat=(-50, 50, 2)).squeeze()

   plt.plot(gpi_zm.latitudes, gpi_zm)
   plt.show()


def pi_xr_ux(*piargs):
   return xr.apply_ufunc(tcpi, *piargs,\
            input_core_dims=[[], [], ['plev'], ['plev'], ['plev']],\
            output_core_dims=[[] for _ in range(len(piargs))],\
            vectorize=True, dask='parallelized',\
            output_dtypes=[float, float, int, float, float])

if __name__ == '__main__':
   main()
