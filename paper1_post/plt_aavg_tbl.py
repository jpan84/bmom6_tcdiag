import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

DIRIS = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250702_unseed2hPa6m/atm',
         '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm',
         '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm',
         '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/atm',
         '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm']
ALIS = ['UNSEED_90', 'UNSEED_50', 'CTL', 'SEED_50', 'SEED_150']
CTLIX = 2
FILI = 'AAVG_5.0_35.0.nc'

SELMO = np.arange(6, 12)
MOWGTS = np.array([30, 31, 31, 30, 31, 30], dtype=np.int_)

PLTVARS = ['DSE', 'LE', 'KE', 'AM', 'TAUAM', 'PS', 'TS', 'PRECT', 'QFLX', 'SHU', 'AHSRC', 'FSNT', 'FLNT', 'SWCF', 'LWCF', 'FSNS', 'FLNS', 'CLD_BL', 'CLD_FT']
ALTNMS = dict(AHSRC='$Q_a$', QFLX='E')
MULTBY = dict(PRECT=1e3 * 86400, QFLX=86400, PS=1e-2)

def main():
   dss = [xr.open_dataset(os.path.join(dr, FILI)) for dr in DIRIS]
   #print(dss[0].time.dt.month)

   tsel = [ds.sel(time=ds['time'].dt.month.isin(SELMO)) for ds in dss]
   momean = [ds.groupby('time.month').mean() for ds in tsel]
   finvals = [(ds * MOWGTS).sum() / MOWGTS.sum() for ds in momean]

   print(finvals[0])

if __name__ == '__main__':
   main()
