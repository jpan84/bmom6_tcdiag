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

ZMLATS = (-90, 90, 0.5)
LATLAB = np.arange(-90, 90.1, 30)

def main():
   #TODO: AHU integrates to 0
   #TODO: AHT from no-storage assumption
   #TODO: OHU_CAM, OHU_MOM, and hfds are equal
   #TODO: OHU integrates to RESTOM
   #TODO: atmo ener bud approx from LW and PREC

   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   camds = ux.open_mfdataset(camgrid, os.path.join(ARCHV, CASE, HIST_CAM))
   #momds = xr.open_mfdataset(os.path.join(ARCHV, CASE, HIST_MOM))

   print(camds['FSNT'].time.dt)
   ahu = derive_flx(camds, AHU)
   print(ahu.time.dt)
   #exit()
   ahu_int = wtm(cam_glob_integ(ahu, camds.uxgrid.face_areas))
   plot_cam_zm(ahu, 'AHU: %.2f W m$^{-2}$' % ahu_int, 'ahu.png')

   restom = derive_flx(camds, RESTOM)
   restom_int = wtm(cam_glob_integ(restom, camds.uxgrid.face_areas))
   plot_cam_zm(restom, 'RESTOM: %.2f W m$^{-2}$' % restom_int, 'restom.png')

def derive_flx(ds, varsfacs):
   return sum([vf[0] * ds[vf[1]] for vf in varsfacs])

def cam_glob_integ(da, areas):
   return (da * areas).sum(dim='n_face') / areas.sum(dim='n_face')

def plot_cam_zm(da, title, outfil):
   toplt = da.zonal_mean(lat=ZMLATS)
   toplt = toplt.assign_coords(coords=dict(time=da.time))
   #print(toplt.time)
   #print(toplt.time.dt)
   sinlat = np.sin(np.deg2rad(toplt.latitudes))
   print(toplt.shape)
   print(wtm(toplt).shape)
   plt.plot(sinlat, wtm(toplt))
   plt.xlabel('Lat [°]')
   plt.title(title)
   plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
   plt.savefig(os.path.join(DIRO, outfil))
   plt.close()

if __name__ == '__main__':
   main()
