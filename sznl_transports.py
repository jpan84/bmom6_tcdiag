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
lncolors = plt.cm.jet(np.linspace(0, 1, 12))
monmeans = lambda da: da.groupby('time.month').mean()

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   print('Opening history files...')
   momdss = [xr.open_mfdataset(os.path.join(ARCHV, CASE1, MOMH))]
   #camdss = [xr.open_mfdataset(os.path.join(ARCHV, CASE1, CAMH))]
   if CASE2 is None:
      momdss.append(xr.zeros_like(momdss[0]))
      #camdss.append(xr.zeros_like(camdss[0]))
   else:
      momdss.append(xr.open_mfdataset(os.path.join(ARCHV, CASE2, MOMH)))
      #camdss.append(xr.open_mfdataset(os.path.join(ARCHV, CASE2, CAMH)))
   #camdss = [ds.rename(dict(lat='lat_h')) for ds in camdss]
   momdss = [ds.rename(dict(yh='lat_h')) for ds in momdss]

   print('Computing OHT...')
   oht = [monmeans((ds.T_ady_2d + ds.T_diffy_2d).sum(dim='xh')) for ds in momdss]
   oht = oht[1] - oht[0]
   #plt.plot(oht.yq, oht.mean(dim='time').values)
   #plt.show()

   sinlat = np.sin(np.deg2rad(oht.yq))
   for mo in oht.month:
      plt.plot(sinlat, oht.sel(month=mo), label=mo.values, color=lncolors[mo.values-1])
   plt.xlabel('Lat [°]')
   plt.ylabel('OHT')
   plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=8))
   plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB)
   plt.savefig(os.path.join(DIRO, '%s.oht_mon.png' % CASE1.split('.')[-1], bbox_inches='tight'))
   plt.show()
   plt.close()

if __name__ == '__main__':
   main()
