#Joshua Pan orig from Oct 2023 for docn runs
#plot and compare zonal-mean time-averaged quantities

#TOPDIR = '/glade/scratch/zarzycki/archive/'
#DIR1 = 'QPC5-ne30np4-aquap10-seed3x3/atm/hist'
#DIR2 = 'QPC5-ne30np4-aquap10-unseed/atm/hist'
#FN = 'QPC5-ne30np4-aquap10-*.cam.h0.*regrid.nc'
OUTDIR = 'linevslat/'
HISTDIMS = set(['time', 'lat', 'lon'])

import proj3
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pltsettings

#print(dir(proj3))
ARCHV = proj3.ARCHV
CASE = proj3.CASE
H0 = proj3.H0
HISTS = proj3.HISTS
tdel = proj3.tdel
wtm = proj3.weighted_temporal_mean

LATLAB = np.array([-90., -60., -30., 0., 30., 60., 90.])

def main():
   pltsettings.set1()

   print('Opening datasets...')
   ds1 = xr.open_mfdataset(os.path.join(ARCHV, CASE, HISTS, H0))
   ds1 = ds1.assign_coords(coords=dict(time=ds1.time - tdel(days=1))) #timestamp of e.g., 0001-01.nc is Feb 1, so subtract a month from time coord
   ds1 = ds1.assign(dict(U200=ds1.U.sel(lev=200, method='nearest')))

   for dv in ds1.data_vars:
      if set(ds1[dv].dims) == HISTDIMS:
         print('Plotting %s...' % dv)
         annmeans = wtm(ds1, dv)
         line1 = annmeans.mean(dim=['time', 'lon'])
         sinlat = np.sin(np.deg2rad(line1.lat))
         plt.plot(sinlat, line1.values)
         plt.xlabel('Lat [°]')
         plt.ylabel(dv)
         plt.gca().set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB)
         plt.savefig(os.path.join(OUTDIR, '%s.png' % dv), bbox_inches='tight')
         plt.close()

         '''
         ax0 = plt.gca().twinx()
         ax0.plot(sinlat, (line2-line1).values, color='green')
         ax0.axhline(y=0, color='black', linestyle='--')
         ax0.set_ylabel('Difference (unseed-seed)')
         plt.savefig(os.path.join(OUTDIR, '%s.png' % dv), bbox_inches='tight')
         plt.close()
         '''

   print('linevslat.py done.')

if __name__ == '__main__':
   main()
