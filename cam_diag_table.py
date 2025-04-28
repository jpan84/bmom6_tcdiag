import glob
import os
import xarray as xr
import pandas as pd

ARCHV = '/glade/derecho/scratch/jpan/archive/%s/atm/hist'
CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed'
tapes = ['h0a', 'h0i', 'h1i', 'h2a', 'h2i']

master = None #pandas dataframe
idx = 0

def main():
   global master, idx
   for tp in tapes:
      fili = glob.glob(os.path.join(ARCHV % CASE, '*.%s.*.nc' % tp))[0]
      ds = xr.open_dataset(fili)

      for dv in ds.data_vars:
         print('\n', dv)
         dct = ds[dv].attrs
         dct['var'] = dv
         dct['dims'] = ' '.join(ds[dv].dims)
         dct['source'] = ds.source
         dct['htape'] = tp
         dct['time_period_freq'] = ds.attrs['time_period_freq']
         #print(dct)

         if master is None:
            master = pd.DataFrame(dct, index=[idx])
         else:
            master = pd.concat([master, pd.DataFrame(dct, index=[idx])])
         idx += 1

   master.to_csv('cam_diag_vars_%s.csv' % CASE)

if __name__ == '__main__':
   main()
