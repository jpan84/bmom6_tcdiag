import os
import numpy as np
import pandas as pd
import consts as c

DIRI = '~/aquaptc/bmom6_tcdiag/seed_stats'
EVNTS = ['250702_unseed_2hPa6m_unseed_events.parquet', '250415_unseed_production_unseed_events.parquet', '251229_seed_match_seed_events.parquet', '250416_seed1x1_production_seed_events.parquet']

def main():
   dfs = [pd.read_parquet(os.path.join(DIRI, ev)) for ev in EVNTS]
   for ii, df in enumerate(dfs):
      if 'dp [hPa]' in df.columns:
         dfs[ii] = df.rename(columns={'dp [hPa]': 'dp'})
      #   dfs[ii].dp *= 100.
      if 'dp' in df.columns:
         dfs[ii].dp /= 100.
      dfs[ii].rp /= 1e3
      dfs[ii]['clat_abs'] = np.abs(dfs[ii]['clat'])


   print([df.head() for df in dfs])
   dp_stats = [(df['dp'].quantile(.05), df['dp'].mean(), df['dp'].quantile(.95)) for df in dfs]
   dp_fstrs = [rf'{tup[1]:.2f} ({tup[0]:.2f}, {tup[2]:.2f})' for tup in dp_stats]
   rp_stats = [(df['rp'].quantile(.05), df['rp'].mean(), df['rp'].quantile(.95)) for df in dfs]
   rp_fstrs = [rf'{tup[1]:.1f} ({tup[0]:.1f}, {tup[2]:.1f})' for tup in rp_stats]
   clat_stats = [(df['clat_abs'].quantile(.05), df['clat_abs'].mean(), df['clat_abs'].quantile(.95)) for df in dfs]
   clat_fstrs = [rf'{tup[1]:.1f} ({tup[0]:.1f}, {tup[2]:.1f})' for tup in clat_stats]

   print(dp_fstrs)
   print(rp_fstrs)

   #ltx_cols = [{'$dp$ [hPa]': dp_fstrs[ii], '$r_p$ [km]': rp_fstrs[ii] for ii in range(len(EVNTS))}
   #ltx_df = pd.DataFrame(columns=
   #df['$dp$ [hPa]'] = []

   data = {
        c.ALI_LTX[ii]: [dp_fstrs[ii], rp_fstrs[ii], clat_fstrs[ii]] 
        for ii in range(len(EVNTS))
   }

   # 29. Create the DataFrame and set the parameter index
   ltx_df = pd.DataFrame(data, index=[r'$dp$ [hPa]', r'$r_p$ [km]', r'$\phi_c$ [°]'])

   # 30. Insert the empty CTL column at index 2
   ltx_df.insert(2, 'CTL', '---')

   distros = [
        'Natural (DetectNodes)', 
        'Natural (DetectNodes)', 
        '---', 
        'Matched to UNSEED\_50', 
        r'Lat: $\mathcal{U}(5°, 20°)$, $d_p: \mathcal{U}(15 hPa, 40 hPa)$, RMW: \mathcal{U}(150 km, 450 km)$'
    ]
    
    # Use .loc to add the row at a specific index name
   ltx_df.loc['Parameter sampling distributions'] = distros

   # Display and Export
   print(ltx_df)
   latex_table = ltx_df.to_latex(
        escape=False, 
        index=True, 
        column_format='l|ccccc',
        caption="Vortex perturbation parameters for each experiment.",
        label="tab:vortex_params"
   )
   print(latex_table)

if __name__ == '__main__':
   main()
