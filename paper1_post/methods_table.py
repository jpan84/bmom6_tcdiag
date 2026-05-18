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
   dp_fstrs = [rf'\makecell{{{tup[1]:.2f} \\ ({tup[0]:.2f}, {tup[2]:.2f})}}' for tup in dp_stats]
   rp_stats = [(df['rp'].quantile(.05), df['rp'].mean(), df['rp'].quantile(.95)) for df in dfs]
   rp_fstrs = [rf'\makecell{{{tup[1]:.1f} \\ ({tup[0]:.1f}, {tup[2]:.1f})}}' for tup in rp_stats]
   clat_stats = [(df['clat_abs'].quantile(.05), df['clat_abs'].mean(), df['clat_abs'].quantile(.95)) for df in dfs]
   clat_fstrs = [rf'\makecell{{{tup[1]:.1f} \\ ({tup[0]:.1f}, {tup[2]:.1f})}}' for tup in clat_stats]

   print(dp_fstrs)
   print(rp_fstrs)

   intcnt = [df[~df['dp'].isna()].shape[0] / ((df['dt'].iloc[-1] - df['dt'].iloc[0]).days / 365.) for df in dfs]
   #intcnt.insert(2, 0.)
   int_fstrs = [rf'\makecell{{{ic:.1f}}}' for ic in intcnt]

   #ltx_cols = [{'$dp$ [hPa]': dp_fstrs[ii], '$r_p$ [km]': rp_fstrs[ii] for ii in range(len(EVNTS))}
   #ltx_df = pd.DataFrame(columns=
   #df['$dp$ [hPa]'] = []

   data = {
        c.ALI_LTX[ii]: [dp_fstrs[ii], rp_fstrs[ii], clat_fstrs[ii], int_fstrs[ii]] 
        for ii in range(len(EVNTS))
   }

   # 29. Create the DataFrame and set the parameter index
   ltx_df = pd.DataFrame(data, index=[r'$dp$ [hPa]', r'$r_p$ [km]', r'$|\phi_c|$ [°]', 'Annual intervention count'])

   # 30. Insert the empty CTL column at index 2
   ltx_df.insert(2, 'CTL', '---')

   distros = [
        r'\makecell{{Natural\\(DetectNodes)}}', 
        r'\makecell{{Natural\\(DetectNodes)}}', 
        '---', 
        r'\makecell{{Matched to \\UNSEED\_50}}', 
        r'\makecell{{Lat: $\mathcal{U}(5°, 20°)$,\\$dp: \mathcal{U}(15 hPa, 40 hPa)$,\\RMW: $\mathcal{U}(150 km, 450 km)$}}'
    ]
   descr = [
        r'\makecell{{Unseed all TCs\\in warm season\\DetectNodes thresholds\\SLP: 4 hPa → 2 hPa,\\DZ: 15 m → 6 m}}', 
        r'\makecell{{Unseed all TCs\\in warm season\\Default online\\DetectNodes thresholds\\plus $\zeta$ threshold:\\$8 \times 10^{-5}$ s$^{-1}$}}', 
        r'\makecell{{Free-running}}', 
        r'\makecell{{Annual seed count\\and vortex parameters\\matched to \\UNSEED\_50}}', 
        r'\makecell{{1 seed per day\\in warm season}}'
    ]
   intint = ['---' if ii == 2 else '24 hours' for ii in range(5)]
   #print(dfs[0]['dt'])

   ##TODO: add q-factor/moisture treatment
   ##TODO: add actual description (i.e., unseed all TCs every restart)
   ##TODO: add intervention counts (# of un/seeds per year)    
    # Use .loc to add the row at a specific index name
   ltx_df.loc['Parameter sampling distributions'] = distros
   ltx_df.loc['Description'] = descr
   ltx_df.loc['Intervention interval'] = intint
   ltx_df.loc[r'$q$ factor'] = [2.5, 2.5, '---', 2.5, 0.]
   #ltx_df.loc['Annual intervention count'] = intcnt

   new_order = ['Description', 'Intervention interval', 'Annual intervention count', 'Parameter sampling distributions', r'$dp$ [hPa]', r'$r_p$ [km]', r'$|\phi_c|$ [°]', r'$q$ factor']
   ltx_df = ltx_df.reindex(new_order)

   # Display and Export
   print(ltx_df, '\n', '\n')
   latex_table = ltx_df.to_latex(
        escape=False, 
        index=True, 
        column_format='p{2cm}|XXXXX',
        caption="Vortex perturbation parameters for each experiment.",
        label="tab:vortex_params"
   )
   latex_table = latex_table.replace(r'\begin{tabular}', r'\begin{tabularx}{\textwidth}')
   latex_table = latex_table.replace(r'\end{tabular}', r'\end{tabularx}')

   print(latex_table, '\n\n')

   final_output = (
     r"\renewcommand{\arraystretch}{2.8}" + "\n"
     + latex_table
   )
   print(final_output)

if __name__ == '__main__':
   main()
