#plot oht and explore variability
#import warnings
#warnings.filterwarnings('ignore')

import os
import numpy as np
import xarray as xr # requires >= 0.15.1
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt

ARCHV = '/glade/derecho/scratch/jpan/archive/'
CASE = 'b.e23.BMOM.f09_sx0.66av1.aqua.production.0530ywbranch'
HISTS = 'ocn/hist/'
H0 = r'*mom6.hm_[0-9]*.nc'

LATLAB = np.array([-90., -60., -30., 0., 30., 60., 90.])
lncolors = plt.cm.jet(np.linspace(0, 1, 12))

def main():
   print('Opening history files...')
   ds = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H0)).sum(dim='xh')

   oht = ds.T_ady_2d + ds.T_diffy_2d
   #plt.plot(oht.yq, oht.mean(dim='time').values)
   #plt.show()

   monmeans = oht.groupby('time.month').mean()
   sinlat = np.sin(np.deg2rad(monmeans.yq))
   for mo in monmeans.month:
      plt.plot(sinlat, monmeans.sel(month=mo), label=mo.values, color=lncolors[mo.values-1])
   plt.xlabel('Lat [°]')
   plt.ylabel('OHT')
   plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, prop=dict(size=8))
   plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB)
   plt.savefig('%s.oht_mon.png' % CASE.split('.')[-1], bbox_inches='tight')
   plt.show()
   plt.close()

if __name__ == '__main__':
   main()
