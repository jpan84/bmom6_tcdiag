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

#FILI = 'AAVG_-90.0_90.0.nc'
#SELMO = np.arange(1, 13)
#MOWGTS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.int_)

PLTVARS = ['DSE', 'LE', 'KE', 'KE_MEAN', 'AM', 'TAUAM', 'PS', 'TS', 'PRECT', 'QFLX', 'SHU', 'AHSRC', 'FSNT', 'FLNT', 'SWCF', 'LWCF', 'FSNS', 'FLNS', 'CLD_BL', 'CLD_FT']
ALTNMS = dict(AHSRC='$Q_a$', QFLX='E')
MULTBY = dict(PRECT=1e3 * 86400, QFLX=86400, PS=1e-2)

VARLBLS = dict(TS='SST [K]', DSE='DSE [J m$^{-2}$]', LE='LE [J m$^{-2}$]', PRECT='P [mm d$^{-1}$]', QFLX='E [mm d$^{-1}$]')

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

   #plot containing all vars
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
   plt.close()
   #plt.show()

   #thermo plot to incl in paper 1
   # 1. Define row blocks: (Variable List, Shared Colorbar Label)
   blocks = [
       (['TS'], 'Anom [K]'),
       (['EKE'], 'Anom [J m$^{-2}$]'),
       (['DSE', 'LE'], 'Anom [J m$^{-2}$]'),
       (['PRECT', 'QFLX'], 'Anom [mm d$^{-1}$]')
   ]
   
   # Extract row slices and compute custom EKE mapping (KE anomaly relative to CTL)
   v_row = {pv: valtbl[PLTVARS.index(pv), :] for pv in PLTVARS if pv != 'KE'}
   d_row = {pv: diftbl[PLTVARS.index(pv), :] for pv in PLTVARS if pv != 'KE'}
   v_row['EKE'] = valtbl[PLTVARS.index('KE'), :] - valtbl[PLTVARS.index('KE_MEAN'), :]
   d_row['EKE'] = v_row['EKE'] - v_row['EKE'][CTLIX]

   # 2. Build Grid Topology dynamically with 0.4 height spacers between blocks
   grid_vars = [v for b in blocks for v in b[0] + ['SPACER']][:-1]
   hratios = [0.4 if v == 'SPACER' else 0.7 for v in grid_vars]
   fig, axs = plt.subplots(len(grid_vars), 1, figsize=(8, 2. * len(blocks)), 
                           gridspec_kw={'height_ratios': hratios, 'hspace': 0}) # hspace=0 removes gaps inside blocks

   # 3. Process each group block cleanly
   for vars_in_block, cbar_lbl in blocks:
      block_axs = [axs[grid_vars.index(v)] for v in vars_in_block]
      b_max = max(np.abs(d_row[v]).max() for v in vars_in_block)
      
      for v in vars_in_block:
         ax = axs[grid_vars.index(v)]
         im = ax.imshow(d_row[v][np.newaxis, :], cmap='bwr', aspect='auto', vmin=-b_max, vmax=b_max)
         
         # Map clean y-labels and inject cell text matrix
         ax.set_yticks([0]), ax.set_yticklabels([VARLBLS.get(v, v).split(' ')[0]])
         ax.set_xticks(range(len(ALIS))), ax.set_xticklabels(ALIS if v == vars_in_block[-1] else [])
         for cc, (val, dif) in enumerate(zip(v_row[v], d_row[v])):
            ax.text(cc, 0, ('%+.2e' % dif) if cc != CTLIX else ('%.2e' % val), ha='center', va='center', 
                     color='white' if abs(dif/b_max) > 0.5 else 'black')
      
      fig.colorbar(im, ax=block_axs, orientation='vertical', pad=0.03, label=cbar_lbl)

   # 4. Hide structural spacer rows from the canvas and save
   [ax.axis('off') for name, ax in zip(grid_vars, axs) if name == 'SPACER']
   plt.savefig(FILI + '_clean_blocks.png', dpi=300, bbox_inches='tight'), plt.close()

# --- Radiation and Cloud Plot for Paper 1 ---
   # 1. Define row blocks: (Variable List, Shared Colorbar Label)
   blocks = [
       (['FSNT', 'SWCF', 'FLNT', 'LWCF'], 'Anom [W m$^{-2}$]'),
       (['CLD_BL', 'CLD_FT'], 'Anom [g kg$^{-1}$]')
   ]
   
   # Extract row slices for the targeted variables
   v_row = {pv: valtbl[PLTVARS.index(pv), :] for pv in PLTVARS}
   d_row = {pv: diftbl[PLTVARS.index(pv), :] for pv in PLTVARS}

   # 2. Build Grid Topology dynamically with 0.4 height spacers between blocks and 0.7 data rows
   grid_vars = [v for b in blocks for v in b[0] + ['SPACER']][:-1]
   hratios = [0.4 if v == 'SPACER' else 0.7 for v in grid_vars]
   fig, axs = plt.subplots(len(grid_vars), 1, figsize=(8, 2. * len(blocks)), 
                           gridspec_kw={'height_ratios': hratios, 'hspace': 0}) # hspace=0 removes gaps inside blocks

   # 3. Process each group block cleanly
   for vars_in_block, cbar_lbl in blocks:
      block_axs = [axs[grid_vars.index(v)] for v in vars_in_block]
      b_max = max(np.abs(d_row[v]).max() for v in vars_in_block)
      
      for v in vars_in_block:
         ax = axs[grid_vars.index(v)]
         im = ax.imshow(d_row[v][np.newaxis, :], cmap='bwr', aspect='auto', vmin=-b_max, vmax=b_max)
         
         # Map clean y-labels and inject cell text matrix
         ax.set_yticks([0]), ax.set_yticklabels([VARLBLS.get(v, v).split(' ')[0]])
         ax.set_xticks(range(len(ALIS))), ax.set_xticklabels(ALIS if v == vars_in_block[-1] else [])
         for cc, (val, dif) in enumerate(zip(v_row[v], d_row[v])):
            ax.text(cc, 0, ('%+.2e' % dif) if cc != CTLIX else ('%.2e' % val), ha='center', va='center', 
                     color='white' if abs(dif/b_max) > 0.5 else 'black')
      
      fig.colorbar(im, ax=block_axs, orientation='vertical', pad=0.03, label=cbar_lbl)

   # 4. Hide structural spacer rows from the canvas and save
   [ax.axis('off') for name, ax in zip(grid_vars, axs) if name == 'SPACER']
   plt.savefig(FILI + '_radiation_clouds.png', dpi=300, bbox_inches='tight'), plt.close()

if __name__ == '__main__':
   main()
