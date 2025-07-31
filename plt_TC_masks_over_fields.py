import uxarray as ux
import numpy as np
import pandas as pd
import os
import cftime

CASE = 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl'
ATM = '/glade/derecho/scratch/jpan/archive/%s/atm'
H1I = os.path.join(ATM % CASE, 'hist/') #h1i path
MASKS = os.path.join(ATM % CASE, 'nff_2mps/') #masks path
PARQ = '/glade/u/home/jpan/aquaptc/tempest/250417_ctrl.parquet' #TC parquet

STRF = lambda dtobj: f'*{dtobj.year:04}-{dtobj.month:02}-{dtobj.day:02}-{3600*dtobj.hour:05}*'

def main():
   tcdf = pd.read_parquet(PARQ)
   print(tcdf)

   for ii, row in tcdf.iterrows(): #For timestep in TC parquet

      if np.isnan(row['year']) or row['year'] != 2010.:
         continue
      truedt = cftime.DatetimeNoLeap(int(row['year']) - 2000, int(row['month']), int(row['day']), hour=int(row['hour']))
      print(truedt, STRF(truedt))

      #open h1i
      #open mask

      #For field in ['U850', 'V850', 'PRECT', 'OMEGA500', 'PS']
         #contourf/rasterize plot field
         #contour mask

if __name__ == '__main__':
   main()
