import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FILI = '~/aquaptc/tempest/260415_density_5exp/sznl_climo.csv' #from track_dens.py
COLS = ['mo_in_szn', 'varnm', 'case', 'NH', 'SH', 'totyrs']
nh_warm_mo = [6, 9]
sh_warm_mo = [2, 3]

DIRE = '~/aquaptc/tempest/'
EVNTS = [('250702_unseed_2hPa6m_unseed_events.parquet', 'us'), ('250415_unseed_production_unseed_events.parquet', 'us'), None, ('251229_seed_match_seed_events.parquet', 'sd'), ('250416_seed1x1_production_seed_events.parquet', 'sd')]
ORDER = ['unseed2', 'unseed', 'ctrl', 'mseed', 'seed']

def main():
   climo = pd.read_csv(FILI, header=None, names=COLS).drop_duplicates()
   climo[['NH', 'SH']] = climo[['NH', 'SH']].div(climo['totyrs'], axis=0)
   print(climo)

   evdfs = [pd.read_parquet(os.path.join(DIRE, ev[0])) if ev is not None else None for ev in EVNTS]
   #TODO: filter unseed dfs by (~dp.isnan())
   evdfs = [df[~df['rp'].isna()] if df is not None else None for df in evdfs]
   print([df.shape if df is not None else None for df in evdfs])

   #print(tot_nTCs(climo))
   for ii, df in enumerate(evdfs):
      print('\nWorking on case', ORDER[ii])
      if EVNTS[ii] is None:
         continue
      selcase = climo[climo['case'] == ORDER[ii]]
      ev_pyr = df.shape[0] / selcase['totyrs'].iloc[0] #un/seed events per year
      tccnt = selcase[selcase['varnm'] == 'genesis points (storm count)']
      nTCs = (tccnt['NH'] + tccnt['SH']).sum()
      dn = nTCs - tot_nTCs(climo[climo['case'] == 'ctrl']).iloc[0]
      #print(dn)
      if EVNTS[ii][1] == 'sd':
         seed_gen = selcase[selcase['varnm'] == 'seeded genesis points']
         n_converted = (seed_gen['NH'] + seed_gen['SH']).sum()
         print('The seed conversion rate for %s is %.3f' % (ORDER[ii], n_converted / ev_pyr))
         print('The efficacy of seeds for %s is %.3f' % (ORDER[ii], dn / n_converted))

   exit()


   nh_warm = df[df['mo_in_szn'].isin(nh_warm_mo)].groupby(['varnm', 'case'], sort=False)['NH'].sum()
   sh_warm = df[df['mo_in_szn'].isin(sh_warm_mo)].groupby(['varnm', 'case'], sort=False)['SH'].sum()

   print(nh_warm + sh_warm)

   # ... (Your existing loading and calculation code) ...
   # Assuming combined is the result of nh_warm + sh_warm
   combined = nh_warm + sh_warm
   
   # 1. Filter for just the 3 variables you requested
   target_vars = ['ACE [$10^4$ kt$^2$] ', 'genesis points (storm count)', 'genesis points of hurricanes']
   # Use .loc with slice(None) to select specific index levels
   plot_series = combined.loc[target_vars].rename(index={'genesis points (storm count)': 'Number of TCs', 'genesis points of hurricanes': 'Number of hurricanes'})
   
   # 2. Reshape from Series to Matrix (Rows=Vars, Cols=Cases)
   # unstack('case') moves the 'case' index level to columns
   plot_df = plot_series.unstack('case')
   
   # Reorder columns to ensure 'ctrl' is in the middle for visual balance
   case_order = ['unseed2', 'unseed', 'ctrl', 'mseed', 'seed']
   plot_df = plot_df[case_order]

   # 3. Calculate Percent Difference relative to 'ctrl'
   pct_diff = plot_df.div(plot_df['ctrl'], axis=0).subtract(1).mul(100)

   # 4. Plotting
   fig, ax = plt.subplots(figsize=(12, 7))
   
   # Center cmap at 0 using vmin/vmax
   # Setting vmin/vmax to +/- 200% for clear color gradients
   im = ax.imshow(pct_diff, cmap='bwr', aspect='auto', vmin=-200, vmax=200)

   # Add text annotations
   for i in range(len(plot_df.index)):
      for j in range(len(plot_df.columns)):
         val = plot_df.iloc[i, j]
         diff = pct_diff.iloc[i, j]
           
         if plot_df.columns[j] == 'ctrl':
            label = f"{val:.2f}\n(REF)"
         else:
            prefix = "+" if diff > 0 else ""
            label = f"{val:.2f}\n({prefix}{diff:.1f}%)"
           
         # Change text color based on background darkness
         text_col = "white" if abs(diff) > 130 else "black"
         ax.text(j, i, label, ha="center", va="center", color=text_col, fontweight='bold')

   # Styling
   ax.set_xticks(np.arange(len(case_order)))
   ax.set_xticklabels(case_order)
   ax.set_yticks(np.arange(len(plot_df.index)))
   ax.set_yticklabels(plot_df.index)
   
   plt.title("Warm Season Climatology: Absolute Values and % Change vs Ctrl", pad=20)
   plt.colorbar(im, label="Percent Difference vs Ctrl (%)")
   plt.tight_layout()
   plt.show()

#get total number of TCs by case
def tot_nTCs(df):
   tccnt = df[df['varnm'] == 'genesis points (storm count)']
   return (tccnt['NH'] + tccnt['SH']).groupby(df['case'], sort=False).sum()

if __name__ == '__main__':
   main()
