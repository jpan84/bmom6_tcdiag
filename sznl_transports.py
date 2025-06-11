import os
import numpy as np
import xarray as xr
#from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt

import sznl_zm_ux

wmat = sznl_zm_ux.wmat
print(wmat)

DIRO = 'sznl_transports_250417_ctrl'
ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE1 = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl'
CASE2 = None
MOMH = 'ocn/hist/*mom6.hm*[0-9][0-9][0-9][0-9]-[0-9][0-9].nc'
CAMH = 'atm/hist_regrid_mom_onpres/*h0a*.nc'

LATLAB = np.array([-90., -60., -30., 0., 30., 60., 90.])
SZNS = sznl_zm_ux.SZNS
lncolors = plt.cm.jet(np.linspace(0, 1, 4))
monmeans = lambda da: da.groupby('time.month').mean()
monzm = lambda da: monmeans(da.mean(dim='lon'))

g = 9.81
cp = 1004
lv = 2500840
a = 6.371e6

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   print('Opening history files...')
   momdss = [xr.open_mfdataset(os.path.join(ARCHV, CASE1, MOMH))]
   camdss = [xr.open_mfdataset(os.path.join(ARCHV, CASE1, CAMH))]
   if CASE2 is None:
      momdss.append(xr.zeros_like(momdss[0]))
      camdss.append(xr.zeros_like(camdss[0]))
   else:
      momdss.append(xr.open_mfdataset(os.path.join(ARCHV, CASE2, MOMH)))
      camdss.append(xr.open_mfdataset(os.path.join(ARCHV, CASE2, CAMH)))
   camdss = [ds.rename(dict(lat='lat_h')) for ds in camdss]
   momdss = [ds.rename(dict(yh='lat_h')) for ds in momdss]

   print('Computing OHT...')
   oht = [monmeans((ds.T_ady_2d + ds.T_diffy_2d).sum(dim='xh')) for ds in momdss]
   oht = oht[0] - oht[1]
   oht_sznl = wmat @ oht.data
   #plt.plot(oht.yq, oht.mean(dim='time').values)
   #plt.show()
   sinlat_q = np.sin(np.deg2rad(oht.yq))
   plt_sznl(sinlat_q, oht_sznl, 'OHT', '[W]')

   print('Computing AHT...')
   gpf = [g * monzm(ds['V']) * monzm(ds['Z3']) for ds in camdss]
   shf = [cp * monzm(ds['VT']) for ds in camdss]
   lhf = [lv * monzm(ds['VQ']) for ds in camdss]
   mseflx = sum([f[0] - f[1] for f in [gpf, shf, lhf]])
   latcirc = 2 * np.pi * a * np.cos(np.deg2rad(mseflx.lat_h))
   aht = -latcirc / g * mseflx.integrate('plev')
   aht = aht.transpose('month', ...)
   print(aht.shape)
   aht_sznl = wmat @ aht.data
   sinlat_cam_h = np.sin(np.deg2rad(mseflx.lat_h))
   plt_sznl(sinlat_cam_h, aht_sznl, 'AHT', '[W]')

   print('Computing AMT...')
   umf = [monzm(ds['UV']) for ds in camdss]
   umf = umf[0] - umf[1]
   amt = -latcirc / g * umf.integrate('plev')
   amt = amt.transpose('month', ...)
   amt_sznl = wmat @ amt.data
   plt_sznl(sinlat_cam_h, amt_sznl, 'AMT', '[N]')

   print('Computing OMT...')
   omf = [monmeans(ds['rhoinsitu'] * ds['uv']).mean(dim='xh') for ds in momdss]
   omf = omf[0] - omf[1]
   omt = latcirc * omf.integrate('zl')
   omt_sznl = wmat @ omt.data
   sinlat_h = np.sin(np.deg2rad(omt.yh))
   plt_sznl(sinlat_h, omt_sznl, 'OMT', '[N]')

   print('Done.')

#vals should have shape (time, lat)
def plt_sznl(sinlats, vals, name, units, linestyle='solid', close=True):
   for tt in range(vals.shape[0]):
      plt.plot(sinlats, vals[tt, :], label=SZNS[tt], color=lncolors[tt], linestyle=linestyle)
   plt.xlabel('Lat [°]')
   plt.ylabel(name + ' ' + units)
   plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=10))
   plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB.astype(np.int_))
   plt.savefig(os.path.join(DIRO, '%s.%s.png' % (CASE1.split('.')[-1], name)), bbox_inches='tight')
   if close:
      plt.close()

if __name__ == '__main__':
   main()
