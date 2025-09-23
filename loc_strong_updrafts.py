import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from dask.diagnostics import ProgressBar

DIRO = './TC_updraft_contrib'
DIRI = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/hist_0012-0014_h1i/'
CASES = ['250415_unseed', '250417_ctrl', '250416_seed1x1']

#OMVALS = 10. ** np.arange(.6, 1.7, .2)
OMVALS = np.arange(16, 25, 2)
TSLC = slice(None, None, 10)
a_e = 6.371e6

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   print('Opening datasets...')
   h1is = [xr.open_mfdataset(os.path.join(DIRI % cs, 'cat_h1i_0012-0014.nc')) for cs in CASES]
   msks = [xr.open_mfdataset(os.path.join(DIRI % cs, 'TC_R4_masked_h1i_0012-0014.nc')) for cs in CASES]

   #subset time for testing
   print('Subsetting time...')
   h1is = [h1.isel(time=TSLC) for h1 in h1is]
   msks = [ms.isel(time=TSLC) for ms in msks]

   '''
   #repeat area along time axis
   print('Repeating area array...')
   arearep = [h1['area'].expand_dims(dim=dict(time=h1.time), axis=0) for h1 in h1is]
   #print([ar.sum(dim='ncol').mean(dim='time').compute() for ar in arearep])

   #idxs where each updraft velocity threshold is met. shape: (case, velocity)
   print('Finding indices with strong updrafts...')
   allidxs = [[np.where(-h1['OMEGA500'] > om) for om in OMVALS] for h1 in h1is]
   tcsidxs = [[np.where(-ms['OMEGA500'] > om) for om in OMVALS] for ms in msks]

   with ProgressBar():
      print('Omega values are', OMVALS)
      totareas = [[ar.isel(time=allidxs[ii][jj][0], ncol=allidxs[ii][jj][1]).sum().compute().data\
                   for jj, om in enumerate(OMVALS)] for ii, ar in enumerate(arearep)]
      tcsareas = [[ar.isel(time=tcsidxs[ii][jj][0], ncol=tcsidxs[ii][jj][1]).sum().compute().data\
                   for jj, om in enumerate(OMVALS)] for ii, ar in enumerate(arearep)]
   '''

   #use binary masks instead of np.where
   allmet = [[(-h1['OMEGA500'] > om) for om in OMVALS] for h1 in h1is]
   tcsmet = [[(-ms['OMEGA500'] > om) for om in OMVALS] for ms in msks]
   #print(allmet[1])

   #
   totareas, tcsareas = None, None
   with ProgressBar():
      print('Omega values are', OMVALS)
      totareas = [[h1['area'].broadcast_like(h1['OMEGA500']).where(met).sum().compute().data\
                     for met in allmet[ii]] for ii, h1 in enumerate(h1is)]
      tcsareas = [[h1['area'].broadcast_like(h1['OMEGA500']).where(met).sum().compute().data\
                     for met in tcsmet[ii]] for ii, h1 in enumerate(h1is)]

   print(totareas[0])

   print('Plotting...')
   ax1 = plt.axes()
   ax2 = ax1.twinx()
   for ii, ta in enumerate(totareas):
      ax1.plot(OMVALS, a_e**2 * np.array(ta) / h1is[ii].time.size)
      ax2.plot(OMVALS, np.array(tcsareas[ii]) / np.array(ta), linestyle='dashed')
   #ax1.set_xscale('log')
   ax1.set_yscale('log')

   plt.savefig(os.path.join(DIRO, 'updraft_area_TC_frac.png'))
   plt.show()

   print('loc_strong_updrafts.py done.')

if __name__ == '__main__':
   main()
