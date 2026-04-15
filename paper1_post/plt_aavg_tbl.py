import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

DIRIS = [
         '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250702_unseed2hPa6m/atm',
         '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm',
         '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm',
         '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/atm',
         '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/atm'
         ]
ALIS = ['UNSEED_90', 'UNSEED_50', 'CTL', 'SEED_50', 'SEED_150']
CTLIX = 2

FILI = 'AAVG_5.0_35.0.nc'
SELMO = np.arange(6, 12)
MOWGTS = np.array([30, 31, 31, 30, 31, 30], dtype=np.int_)

FILI = 'AAVG_-90.0_90.0.nc'
SELMO = np.arange(1, 13)
MOWGTS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.int_)

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

   valtbl = np.zeros((len(PLTVARS), len(ALIS)))
   diftbl = np.zeros_like(valtbl)
   toshade = np.zeros_like(diftbl)
   for rr, pv in enumerate(PLTVARS):
      valtbl[rr, :] = np.array([fv[pv] for fv in finvals]) * MULTBY.get(pv, 1.)
      diftbl[rr, :] = valtbl[rr, :] - valtbl[rr, CTLIX]
      toshade[rr, :] = diftbl[rr, :] / np.abs(diftbl[rr, :]).max()

   fig, ax = plt.subplots(figsize=(6, 18))
   im = ax.imshow(toshade, cmap='bwr', aspect='auto', vmin=-1, vmax=1)

   ax.set_xticks(range(len(ALIS)))
   ax.set_yticks(range(len(PLTVARS)))
   ax.set_xticklabels(ALIS)
   ax.set_yticklabels(PLTVARS)

   for rr in range(valtbl.shape[0]):
      for cc in range(valtbl.shape[1]):
         lbl = '%+.2e' % (diftbl[rr, cc] if cc != CTLIX else valtbl[rr, cc])
         txtclr = 'white' if abs(toshade[rr, cc]) > 0.5 else 'black'
         ax.text(cc, rr, lbl, ha='center', va='center', color=txtclr)

   plt.savefig(FILI + '.png')
   plt.show()

if __name__ == '__main__':
   main()
