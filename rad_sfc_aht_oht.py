import os
import xarray as xr
import uxarray as ux
import numpy as np
import matplotlib.pyplot as plt
from compowake import selarea, latlon_avg
from globavg import weighted_temporal_mean as wtm

DIRO = 'rad_sfc_250417_ctrl'
ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl'
HIST_CAM = 'atm/hist/*.cam.h0a.*.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
HIST_MOM = 'ocn/hist/*mom6.hm*[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'
momgrid = '/glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed.mom6.static.nc'

RESTOM = [(1, 'FSNT'), (-1, 'FLNT')]
AHU = [(1, 'FSNT'), (-1, 'FSNS'), (1, 'FLNS'), (-1, 'FLNT'), (1, 'SHFLX'), (1, 'LHFLX')]
OHU_CAM = [(1, 'FSNS'), (-1, 'FLNS'), (-1, 'SHFLX'), (-1, 'LHFLX')]
OHU_MOM = [(1, 'rsntds'), (1, 'rlntds'), (1, 'hfsso'), (1, 'hflso')]
OHT = [(1, 'T_ady_2d'), (1, 'T_diffy_2d')]
OMU = [(1, 'taux_bot'), (1, 'tauuo')]

ZMLATS = (-90, 90, 0.5)
LATLAB = np.arange(-90, 90.1, 30)
a = 6.371e6

def main():
   #TODO: AHU integrates to 0
   #TODO: AHT from no-storage assumption
   #TODO: OHU_CAM, OHU_MOM, and hfds are equal
   #TODO: OHU integrates to RESTOM
   #TODO: atmo ener bud approx from LW and PREC

   #TODO: OMT from uv
   #TODO: AMT from UV
   #TODO: sum of ocean stresses
   #TODO: sum of atmo stresses

   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   camds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HIST_CAM))
   momds = xr.open_mfdataset(os.path.join(ARCHV, CASE, HIST_MOM))
   geogrid = xr.open_dataset(momgrid)

   ahu = derive_flx(camds, AHU)
   ahu_int = wtm(cam_glob_integ(ahu, camds.uxgrid.face_areas))
   plot_cam_zm(wtm(ahu), 'AHU: %.2f W m$^{-2}$' % ahu_int, 'ahu.png')

   print('RESTOM')
   restom = derive_flx(camds, RESTOM)
   restom_int = wtm(cam_glob_integ(restom, camds.uxgrid.face_areas))
   plot_cam_zm(wtm(restom), 'RESTOM: %.2f W m$^{-2}$' % restom_int, 'restom.png')

   print('OHU CAM')
   ohu_cam = derive_flx(camds, OHU_CAM)
   ohu_cam_int = wtm(cam_glob_integ(ohu_cam, camds.uxgrid.face_areas))
   plot_cam_zm(wtm(ohu_cam), 'OHU: [W m$^{-2}$]', 'ohu_cam.png', label='CAM %.2f' % ohu_cam_int, close=False)

   print('OHU MOM')
   ohu_mom = wtm(derive_flx(momds, OHU_MOM))
   ohu_mom, momlat = selarea(ohu_mom, geogrid, (-90., 90.), (0., 360.)) 
   #print(momlat)
   geolat = momlat.isel(xh=0).data
   ohu_mom_zm = ohu_mom.mean(dim='xh')
   ohu_mom_int = latlon_avg(ohu_mom, momlat)
   plt.plot(np.sin(np.deg2rad(geolat)), ohu_mom_zm, label='MOM %.2f' % ohu_mom_int)
   plt.legend()
   plt.savefig('OHU.png', bbox_inches='tight')
   plt.close()

   print('OHT')
   oht = wtm(derive_flx(momds, OHT)).sum(dim='xh')
   plt.plot(np.sin(np.deg2rad(momlat)), oht)
   plt.savefig('OHT.png', bbox_inches='tight')
   plt.close()

   print('AOHT from RESTOM')
   aoht_restom = hfl2ht(wtm(restom).zonal_mean(lat=geolat))
   plt.plot(aoht_restom.latitudes, aoht_restom)
   plt.title('AOHT from RESTOM: [W]')
   plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
   plt.savefig(os.path.join(DIRO, 'aoht_restom.png'))

   outds = xr.Dataset(data_vars=dict(ahu=ahu, restom=restom, ohu_cam=ohu_cam, ohu_mom=ohu_mom, aoht_restom=aoht_restom))
   outds.to_netcdf(os.path.join(DIRO, 'out.nc'))

   print('Done.')

def derive_flx(ds, varsfacs):
   return sum([vf[0] * ds[vf[1]] for vf in varsfacs])

def cam_glob_integ(da, areas):
   return (da * areas).sum(dim='n_face') / areas.sum(dim='n_face')

def plot_cam_zm(da, title, outfil, label=None, close=True):
   toplt = da.zonal_mean(lat=ZMLATS)
   #toplt = toplt.assign_coords(coords=dict(time=da.time))
   sinlat = np.sin(np.deg2rad(toplt.latitudes))
   plt.plot(sinlat, toplt, label=label)
   plt.xlabel('Lat [°]')
   plt.title(title)
   plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
   plt.savefig(os.path.join(DIRO, outfil))
   if close:
      plt.close()

def cum_int(da, dim, wgts=1):
   intda = xr.zeros_like(da)
   for co in da[dim]:
      #print(co)
      subda = (wgts * da).sel({dim: slice(None, co)})
      intda.loc[{dim: co}] = subda.integrate(coord=dim)
   return intda

#zonal-mean heat flux to meridional heat transport
#uses the zero-storage assumption
def hfl2ht(hfl, lat='latitudes'):
   wgts = np.cos(np.deg2rad(hfl[lat]))
   return 2 * np.pi * a**2 * np.pi / 180 * cum_int(hfl, lat, wgts=wgts)

if __name__ == '__main__':
   main()
