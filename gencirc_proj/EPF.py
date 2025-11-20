import os
import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar

DIRI = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr/atm'
#EHFFIL = '0012_JJASON_onpres_0.25_VTeddy.nc'
#EMFFIL = '0012_JJASON_onpres_0.25_VUeddy.nc'
MMCFIL = '0012_JJASON_onpres_0.25_zm.nc'
TOTFLX = '0012_JJASON_onpres_0.25_zmprods.nc'

Rd = 287
cp = 1005
ae = 6.371e6
grav = 9.81
p0 = 1000 #hPa
Om = 7.29e-5
Hd = 7e3
dens0 = p0 * 100 / grav / Hd

PNM = 'plev'
LATNM = 'lat'

def main():
   inds = xr.open_mfdataset([os.path.join(DIRI, MMCFIL), os.path.join(DIRI, TOTFLX)]).squeeze()
   print(inds)
   if inds[PNM].units == 'Pa':
      inds = inds.assign_coords(coords={PNM: inds[PNM] / 100})
      inds[PNM].attrs['units'] = 'hPa'

   #compute relevant constant scaling factors
   exn = exner_hPa(inds[PNM])
   coslat = np.cos(np.deg2rad(inds['lat']))
   sinlat = np.sin(np.deg2rad(inds['lat']))
   dens = dens0 * inds[PNM] / p0
   EPfac = ae * dens * coslat
   inds = inds.assign_coords(zs=(PNM, -Hd * np.log(inds[PNM].data / p0)))
   print(inds['zs'])

   print('Computing EHF and EMF using Reynolds decomp...')
   vth_tot = inds['VT'] / exn
   th_zm = inds['T'] / exn
   vth_mmc = inds['V'] * th_zm
   ehf = vth_tot - vth_mmc
   emf = inds['VU'] - inds['V'] * inds['U']

   print('Computing vT velocity streamfunc...')
   th_z = th_zm.differentiate('zs', edge_order=2)
   Psi_vT =  ehf / th_z

   print('EPphi, EMF term...')
   epy_emf = EPfac * -emf

   print('EPphi, uw adv...')
   U_z = inds['U'].differentiate('zs', edge_order=2)
   epy_adv = EPfac * Psi_vT * U_z

   print('EPpz, EHF term...')
   epz_ehf = EPfac * Psi_vT * 2 * Om * sinlat

   print('EPpz, uv adv...')
   U_y = yderiv(inds['U'], coslat)
   epz_adv = EPfac * Psi_vT * -U_y

   EPterms = [epy_emf, epy_adv, epz_ehf, epz_adv]

def exner_hPa(plevs):
   assert(plevs.units == 'hPa')
   return (plevs / p0) ** (Rd / cp)

def yderiv(var, coslat):
   return np.rad2deg((var * coslat).differentiate(LATNM) / ae / coslat)

if __name__ == '__main__':
   main()
