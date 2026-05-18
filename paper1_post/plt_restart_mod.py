import os
import consts as c
import pandas as pd
import numpy as np
import uxarray as ux
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter as SclrFmt

DIRO = './unseed_restart_plts'
DIRO = './mseed_restart_plts'

RSDIR = '/glade/derecho/scratch/jpan/unseed_restarts_paper1/'
EVNTS = '~/aquaptc/tempest/250415_unseed_production_unseed_events.parquet'
GRIDF = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
GRIDDELTA = 0.25

DTSTR = '0005-08-07-00000'
CLAT, CLON = 18.00722, 98.79274
RADO = 2. #gcd

DTSTR = '0005-09-03-00000'
CLAT, CLON = 29.45009, 141.5423
RADO = 3.
#
#DTSTR = '0006-07-22-00000'
#CLAT, CLON = 24.18184, 2.792722
#RADO = 4.5
#
#DTSTR = '0006-09-01-00000'
#CLAT, CLON = 17.22682, 143.0426
#RADO = 3.
#
#RSDIR = '/glade/derecho/scratch/jpan/mseed_restarts_paper1/'
#EVNTS = '/glade/u/home/jpan/aquaptc/bmom6_tcdiag/seed_stats/251229_seed_match_seed_events.parquet'
#
#DTSTR = '0005-08-17-00000'
#CLAT, CLON = 19.430189, 290.3278144 - 360
#RADO = 3.
#
#DTSTR = '0006-06-07-00000'
#CLAT, CLON = 23.266556, 187.2309796 - 360
#RADO = 3.
#
#DTSTR = '0006-09-30-00000'
#CLAT, CLON = 20.775221, 22.72787575
#RADO = 5.
#
#DTSTR = '0006-11-14-00000'
#CLAT, CLON = 20.402225, 210.3689465 - 360
#RADO = 3.

RADO = 5. #standardize lon width

MODSTR = '*%s.nc'
ORISTR = '*%s.nc.ORIG.nc'

PLEVS = np.arange(-1e4, 0, 5e2)
PLEVS = np.concatenate((PLEVS, -PLEVS[::-1]))
VLEVS = np.arange(-50, 51, 5)
TLEVS = levels=np.arange(-10, 0, 1)
TLEVS = np.concatenate((TLEVS, -TLEVS[::-1]))
THELEVS = np.arange(330, 381, 5)

def main():
   orids = ux.open_mfdataset(GRIDF, os.path.join(RSDIR, ORISTR % DTSTR)).expand_dims(state=['before'])
   modds = ux.open_mfdataset(GRIDF, os.path.join(RSDIR, MODSTR % DTSTR)).expand_dims(state=['after'])
   modds.uxgrid = orids.uxgrid
   ds = ux.concat([orids, modds], dim='state')
   df = pd.read_parquet(EVNTS)
   dtobj = np.datetime64(DTSTR[:-6])
   myrow = df[df['dt'] == dtobj]
   myrow = df.loc[(df['clat'] - CLAT).abs().idxmin()]
   print(myrow)
   print(ds['lev'].attrs)

   ps = ds['PSDRY'] + ds['dpQ'].sum('lev')
   aiterm = ds['hyai'] * ds['P0']
   biterm = ds['hybi'] * ps
   p_ilev = aiterm + biterm
   dp3d = p_ilev.diff('ilev').rename(dict(ilev='lev')).assign_coords(lev=ds['lev'])
   q = ds['dpQ'] / dp3d
   amterm = ds['hyam'] * ds['P0']
   bmterm = ds['hybm'] * ps
   p_lev = amterm + bmterm
   lev_coord = 1000. * (ds['hyam'] + ds['hybm'])
   ds = ds.assign(q=q, p_lev=p_lev)
   #print(p_ilev / ps)

   CC_bolton = lambda TK: 611.2 * np.exp(17.67 * (TK - 273.15) / (TK - 273.15 + 243.5))
   CC_3410 = lambda TK: 100 * np.exp(-6810.5245 / TK - 5.08984 * np.log(TK) + 55.2966)
   pvap = lambda p, q: p * q / (c.MWH2O / c.MWDRY + q * (1 - c.MWH2O / c.MWDRY))
   rh = 100 * pvap(p_lev, q) / CC_3410(ds['T'])
   #rh = ux.UxDataArray(CC_bolton(ds['T']), uxgrid=ds.uxgrid)
   print("dp3d:", dp3d.min().values, dp3d.max().values)
   print("q:", q.min().values, q.max().values)
   print("p_lev:", p_lev.min().values, p_lev.max().values)
   print("T:", ds['T'].min().values, ds['T'].max().values)

   true_ctr = ds.uxgrid.subset.nearest_neighbor((CLON, CLAT), 1)
   az_mean = lambda da: da.azimuthal_mean((true_ctr.face_lon, true_ctr.face_lat), RADO, GRIDDELTA)

   #test of azimuthal mean p anom and v tangential
   #test_p = mirror_azim_mean(az_mean(p_lev).isel(state=0).squeeze())
   #test_p -= test_p.isel(radius=-1)
   #test_vt = mirror_azim_mean(az_mean(uv_to_tang(ds, 'U', 'V', (CLON, CLAT), radbound=RADO)).isel(state=0).squeeze())
   #print(test_vt)
   ##plt.contour(test_p['radius'], lev_coord.isel(state=0), test_p, levels=np.arange(1.5e4, 9.6e4, 5e3), colors='black')
   #plt.contourf(test_vt['radius'], lev_coord.isel(state=0), test_vt, levels=np.arange(-50, 51, 5), cmap='PRGn')
   #plt.colorbar()
   #plt.contour(test_p['radius'], lev_coord.isel(state=0), test_p, levels=PLEVS, colors='black')
   #plt.ylim(1000, 70)
   #plt.yscale('log')
   #plt.show()

   xsect = lambda da: da.cross_section(start=(CLON - RADO, CLAT), end=(CLON + RADO, CLAT), steps=int(2*RADO/GRIDDELTA))

   ##test of p anom and v wind
   #psec = xsect(p_lev.isel(state=0).squeeze()) - az_mean(p_lev).isel(state=0).squeeze().isel(radius=-1).data[:, None]
   #vsec = xsect(ds['V'].isel(state=0).squeeze())
   #plt.contourf(vsec['lon'], lev_coord.isel(state=0), vsec, levels=np.arange(-50, 51, 5), cmap='PRGn')
   #plt.colorbar()
   #plt.contour(psec['lon'], lev_coord.isel(state=0), psec, levels=PLEVS, colors='black')
   #plt.ylim(1000, 70)
   #plt.yscale('log')
   ##plt.show()
   #plt.close()

   ##test of T anom and moisture fields
   #tsec = xsect(ds['T'].isel(state=0).squeeze()) - az_mean(ds['T']).isel(state=0).squeeze().isel(radius=-1).data[:, None]
   #qsec = xsect(q.isel(state=0).squeeze())
   #rhsec = xsect(rh.isel(state=0).squeeze())
   #thesec = xsect(thetae_bolton(p_lev, ds['T'], q).isel(state=0).squeeze())
   ##plt.contourf(qsec['lon'], lev_coord.isel(state=0), qsec, levels=np.arange(2e-3, 3.1e-2, 2e-3), cmap='YlGnBu')
   ##plt.contourf(rhsec['lon'], lev_coord.isel(state=0), rhsec, cmap='YlGnBu')#, levels=np.arange(5, 101, 5))
   #plt.contourf(thesec['lon'], lev_coord.isel(state=0), thesec, cmap='YlGnBu', levels=np.arange(340, 381, 5))
   #plt.colorbar()
   #plt.contour(tsec['lon'], lev_coord.isel(state=0), tsec, levels=TLEVS, colors='black')
   #plt.ylim(1000, 70)
   #plt.yscale('log')
   ##plt.show()
   #plt.close()

   amb_mean = lambda da: da.azimuthal_mean((true_ctr.face_lon, true_ctr.face_lat), 5, 0.5).isel(radius=-1).squeeze().data[..., None]
   #tsec = xsect(ds['T']).squeeze()
   #print(tsec - amb_mean(ds['T']))

   plt.rcParams['figure.figsize'] = (15, 7)
   subplot_kw = dict(ylim=(1000, 70), yscale='log')
   fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, subplot_kw=subplot_kw)

   psec = xsect(p_lev.squeeze()) - amb_mean(p_lev)
   tsec = xsect(ds['T'].squeeze()) - amb_mean(ds['T'])
   vsec = xsect(ds['V'].squeeze())
   qsec = xsect(ds['q'].squeeze())
   thesec = xsect(thetae_bolton(p_lev, ds['T'], q)).squeeze()



   csf = axes[0][0].contourf(vsec['lon'], lev_coord.isel(state=0), vsec.isel(state=0), levels=VLEVS, cmap='PRGn')
   plt.colorbar(csf)
   axes[0][0].contour(psec['lon'], lev_coord.isel(state=0), psec.isel(state=0), levels=PLEVS, colors='black')
   axes[0][0].set_title('v (shaded, m/s),\np\'(contours, 5 hPa)')


   csf = axes[0][1].contourf(vsec['lon'], lev_coord.isel(state=0), vsec.isel(state=1) - vsec.isel(state=0), levels=VLEVS, cmap='bwr')
   plt.colorbar(csf)
   axes[0][1].contour(psec['lon'], lev_coord.isel(state=0), psec.isel(state=1) - psec.isel(state=0), levels=PLEVS, colors='black')
   axes[0][1].text(CLON, 100, f"dp = {myrow['dp'] / 100.:.1f} hPa\nrp = {myrow['rp'] / 1000.:.1f} km\nzp = {myrow['zp'] / 1000.:.2f} km")

   csf = axes[0][2].contourf(vsec['lon'], lev_coord.isel(state=1), vsec.isel(state=1), levels=VLEVS, cmap='PRGn')
   plt.colorbar(csf)
   axes[0][2].contour(psec['lon'], lev_coord.isel(state=1), psec.isel(state=1), levels=PLEVS, colors='black')

   csf = axes[1][0].contourf(thesec['lon'], lev_coord.isel(state=0), thesec.isel(state=0), levels=THELEVS, cmap='YlGnBu', extend='both')
   plt.colorbar(csf)
   axes[1][0].contour(tsec['lon'], lev_coord.isel(state=0), tsec.isel(state=0), levels=TLEVS, colors='gray')
   axes[1][0].set_title('$\\theta_e$ (shaded, K),\nT\'(contours, 1 K)')

   csf = axes[1][1].contourf(qsec['lon'], lev_coord.isel(state=0), qsec.isel(state=1) - qsec.isel(state=0), levels=np.arange(-2e-2, 2.1e-2, 2e-3), cmap='BrBG')
   plt.colorbar(csf)
   axes[1][1].contour(tsec['lon'], lev_coord.isel(state=0), tsec.isel(state=1) - tsec.isel(state=0), levels=TLEVS, colors='black')
   axes[1][1].set_title('$\delta q$ (shaded, kg/kg)')

   csf = axes[1][2].contourf(thesec['lon'], lev_coord.isel(state=1), thesec.isel(state=1), levels=THELEVS, cmap='YlGnBu', extend='both')
   plt.colorbar(csf)
   axes[1][2].contour(thesec['lon'], lev_coord.isel(state=1), tsec.isel(state=1), levels=TLEVS, colors='gray')

   startchar = 'a'
   startchar = 'g'
   sfmt = SclrFmt()
   sfmt.set_scientific(False)
   sfmt.set_useOffset(False)
   colttls = ['Before', 'Perturbation', 'After']
   for ii, ax in enumerate(axes.ravel()):
      ax.set_title('(' + chr(ord(startchar) + ii) + ')', loc='left')
      ax.yaxis.set_major_formatter(sfmt)
      ax.ticklabel_format(style='plain', axis='y')
      ax.tick_params(labelbottom=True, labelleft=True, right=True, top=True)
      ax.set_xlabel('lon\n$\mathbf{%s}$' % (colttls[ii % 3]))
      ax.set_ylabel(ds['lev'].name)
      ax.set_yticks([100, 200, 300, 500, 700, 1000])

   fig.suptitle(DTSTR + '\n(lat, lon) = %.2f, %.2f' % (CLAT, CLON))
   fig.tight_layout()
   plt.savefig(os.path.join(DIRO, DTSTR + '.svg'))
   plt.show()

def mirror_azim_mean(am):
   flip = am.isel(radius=slice(-1, None, -1))
   flip = flip.assign_coords(radius=-flip['radius'])
   return xr.concat([flip, am], dim='radius')

def uv_to_tang(uxds, unm, vnm, ctrcoord, radbound=None):
   myu, myv, mylat, mylon = uxds[unm], uxds[vnm], uxds['lat'], uxds['lon']
   if radbound is not None:
      bc = lambda da: da.subset.bounding_circle(ctrcoord, radbound + 1, element='face centers')
      myu, myv, mylat, mylon = bc(myu), bc(myv), bc(mylat), bc(mylon)
   clon, clat = ctrcoord
   d1 = np.sin(np.deg2rad(clat)) * np.cos(np.deg2rad(mylat)) - np.cos(np.deg2rad(clat)) * np.sin(np.deg2rad(mylat)) * np.cos(np.deg2rad(mylon) - np.deg2rad(clon))
   d2 = np.cos(np.deg2rad(clat)) * np.sin(np.deg2rad(mylon) - np.deg2rad(clon))
   d = np.maximum(1e-25, np.sqrt(d1 ** 2. + d2 ** 2.))

   ufac = d1 / d
   vfac = d2 / d

   #print(type(ufac))
   #return ux.UxDataArray(ufac**2+vfac**2, uxgrid=mylat.uxgrid)
   #return -myu * ufac
   return myv * vfac + myu * ufac

def thetae_bolton(p, T, q):
   mixr = q / (1 - q)
   pvap = p * q / (c.MWH2O / c.MWDRY + q * (1 - c.MWH2O / c.MWDRY))
   tdew = 243.5 / (17.67 / np.log(pvap / 611.2) - 1) + 273.15
   tlcl = 1 / (1 / (tdew - 56) + np.log(T / tdew) / 800) + 56
   thta_lcl = T * (c.P0 / (p - pvap)) ** (c.kapd) * (T / tlcl) ** (0.28 * mixr)
   return thta_lcl * np.exp((3036 / tlcl - 1.78) * mixr * (1 + 0.448 * mixr))

if __name__ == '__main__':
   main()
