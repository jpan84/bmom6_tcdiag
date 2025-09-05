#plot actual and TC-free SSTA time series
import xarray as xr
import numpy as np
import os
from datetime import timedelta as tdel
import matplotlib.pyplot as plt

MASK = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/dSSTAm0.1_AND_TC_R4_runmax15d.nc'
ORIG = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/yhoureddy_sfc_0012-0014.nc'
NoTC = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/cw_15d.nc'
MASKNM = 'tos'
SSTANM = 'tos'

DIRO = './CTRL_cw_series_15d/'

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   print('Opening...')
   mask = xr.open_mfdataset(MASK)[MASKNM]
   ssta = xr.open_mfdataset(ORIG)[SSTANM]
   notc = xr.open_mfdataset(NoTC)[SSTANM]
   print('Opened.')
   ts, te = notc['time'][0], notc['time'][-1]

   print(mask.shape)

   for ii in range(100):
      while True:
         rndidx = np.random.randint(list(mask.shape))
         #print(mask[*rndidx])
         if mask[*rndidx].data:
            print('Cold wake point found')
            selpt = mask[*rndidx]
            #print(mask[*rndidx])
            ti = max(ts, selpt['time'] + tdel(days=-15))
            tf = min(te, selpt['time'] + tdel(days=15))
            print(ti.values, tf.values)


            ssta_ser = ssta.sel(time=slice(ti, tf), xh=selpt['xh'], yh=selpt['yh'])
            notc_ser = notc.sel(time=slice(ti, tf), xh=selpt['xh'], yh=selpt['yh'])

            plt.plot(ssta_ser, label='SSTA')
            plt.plot(notc_ser, label='noTC')
            plt.legend()
            plt.savefig(os.path.join(DIRO, '%d.png' % ii))
            plt.close()
            break

if __name__ == '__main__':
   main()
