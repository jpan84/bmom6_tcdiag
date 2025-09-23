import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

DIRO = './TC_updraft_contrib'
DIRI = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/hist_0012-0014_h1i/'
CASES = ['250415_unseed', '250417_ctrl', '250416_seed1x1']

OMVALS = 10. ** np.arange(.5, 2.1, .5)
TSLC = slice(0, 10)
a_e = 6.371e6

def main():
   if not os.path.exists(DIRO):
      os.makedirs(DIRO)

   h1is = [xr.open_mfdataset(os.path.join(DIRI % cs, 'cat_h1i_0012-0014.nc')) for cs in CASES]
   msks = [xr.open_mfdataset(os.path.join(DIRI % cs, 'TC_R4_masked_h1i_0012-0014.nc')) for cs in CASES]

   #subset time for testing
   h1is = [h1.isel(time=TSLC) for h1 in h1is]
   msks = [ms.isel(time=TSLC) for ms in msks]

   #repeat area along time axis
   arearep = [h1['area'].expand_dims(dim=dict(time=h1.time), axis=0) for h1 in h1is]

   #idxs where each updraft velocity threshold is met. shape: (case, velocity)
   allidxs = [[np.where(-h1['OMEGA500'] > om) for om in OMVALS] for h1 in h1is]
   tcsidxs = [[np.where(-ms['OMEGA500'] > om) for om in OMVALS] for ms in msks]

   totareas = [[ar.isel(time=allidxs[ii][jj][0], ncol=allidxs[ii][jj][1]).sum().compute().data\
                   for jj, om in enumerate(OMVALS)] for ii, ar in enumerate(arearep)]
   tcsareas = [[ar.isel(time=tcsidxs[ii][jj][0], ncol=tcsidxs[ii][jj][1]).sum().compute().data\
                   for jj, om in enumerate(OMVALS)] for ii, ar in enumerate(arearep)]

   ax1 = plt.axes()
   ax2 = ax1.twinx()
   for ii, ta in enumerate(totareas):
      ax1.plot(OMVALS, ta)
      ax2.plot(OMVALS, tcsareas[ii])
   ax1.set_xscale('log')
   ax1.set_yscale('log')

   plt.savefig(os.path.join(DIRO, 'updraft_area_TC_frac.png'))
   plt.show()

if __name__ == '__main__':
   main()
