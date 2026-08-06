import os
import numpy as np
import pandas as pd
import consts as c

DIRI = './TC_preprocess/seed_stats'
EVNTS = [
    '250702_unseed_2hPa6m_events.parquet', 
    '250415_unseed_production_events.parquet', 
    '251229_seed_match_events.parquet', 
    '250416_seed1x1_production_events.parquet'
]

def main():
    dfs = [pd.read_parquet(os.path.join(DIRI, ev)) for ev in EVNTS]
    for ii, df in enumerate(dfs):
        if 'dp [hPa]' in df.columns:
            dfs[ii] = df.rename(columns={'dp [hPa]': 'dp'})
        if 'dp' in df.columns and df['dp'].max() > 500:
            dfs[ii].dp /= 100.
        dfs[ii].rp /= 1e3
        dfs[ii]['clat_abs'] = np.abs(dfs[ii]['clat'])

    print([df.head() for df in dfs])
    
    # Replaced \makecell with base LaTeX \newline
    dp_stats = [(df['dp'].quantile(.05), df['dp'].mean(), df['dp'].quantile(.95)) for df in dfs]
    dp_fstrs = [rf'{tup[1]:.2f} \newline ({tup[0]:.2f}, {tup[2]:.2f})' for tup in dp_stats]
    
    rp_stats = [(df['rp'].quantile(.05), df['rp'].mean(), df['rp'].quantile(.95)) for df in dfs]
    rp_fstrs = [rf'{tup[1]:.1f} \newline ({tup[0]:.1f}, {tup[2]:.1f})' for tup in rp_stats]
    
    clat_stats = [(df['clat_abs'].quantile(.05), df['clat_abs'].mean(), df['clat_abs'].quantile(.95)) for df in dfs]
    clat_fstrs = [rf'{tup[1]:.1f} \newline ({tup[0]:.1f}, {tup[2]:.1f})' for tup in clat_stats]

    print(dp_fstrs)
    print(rp_fstrs)

    intcnt = [df[~df['dp'].isna()].shape[0] / ((df['dt'].iloc[-1] - df['dt'].iloc[0]).days / 365.) for df in dfs]
    int_fstrs = [rf'{ic:.1f}' for ic in intcnt] # No need for \makecell here
    print('The UNSEED_EX years is', dfs[0][~dfs[0]['dp'].isna()].shape[0], (dfs[0]['dt'].iloc[-1] - dfs[0]['dt'].iloc[0]).days / 365.)

    data = {
        c.ALI_LTX[ii]: [dp_fstrs[ii], rp_fstrs[ii], clat_fstrs[ii], int_fstrs[ii]]
        for ii in range(len(EVNTS))
    }

    # Create the DataFrame and set the parameter index
    ltx_df = pd.DataFrame(data, index=[r'$dp$ [hPa]', r'$r_p$ [km]', r'$|\phi_c|$ [°]', 'Annual intervention count'])

    # Insert the empty CTL column at index 2
    ltx_df.insert(2, 'CTL', '---')

    # Swapped \\ for \newline and removed \makecell
    distros = [
        r'Natural \newline (DetectNodes)',
        r'Natural \newline (DetectNodes)',
        '---',
        r'Matched to \newline UNSEED',
        r'Lat: $\mathcal{U}(5°, 20°)$, \newline $dp: \mathcal{U}(15 \text{ hPa}, 40 \text{ hPa})$, \newline RMW: $\mathcal{U}(150 \text{ km}, 450 \text{ km})$'
    ]
    descr = [
        r'Unseed all TCs \newline in warm season \newline DetectNodes thresholds \newline SLP: 4 hPa $\rightarrow$ 2 hPa, \newline DZ: 15 m $\rightarrow$ 6 m',
        r'Unseed all TCs \newline in warm season \newline Default online \newline DetectNodes thresholds \newline plus $\zeta$ threshold: \newline $8 \times 10^{-5}$ s$^{-1}$',
        r'Free-running',
        r'Annual seed count \newline and vortex parameters \newline matched to \newline UNSEED',
        r'1 seed per day \newline in warm season'
    ]
    purp = [
        r'Aggressive unseeding \newline to minimize \newline TC activity',
        r'Moderate unseeding',
        r'Obtain a baseline \newline climatology',
        r'Mirror the forcing \newline in UNSEED',
        r'Aggressive seeding \newline to probe upper bound'
    ]
    
    intint = ['---' if ii == 2 else '24 hours' for ii in range(5)]

    # Use .loc to add the row at a specific index name
    ltx_df.loc['Parameter sampling distributions'] = distros
    ltx_df.loc['Purpose'] = purp
    ltx_df.loc['Description'] = descr
    ltx_df.loc['Intervention interval'] = intint
    ltx_df.loc[r'$q$ factor'] = ['2.5', '2.5', '---', '2.5', '0']

    new_order = [
        'Purpose', 'Description', 'Intervention interval', 'Annual intervention count', 
        'Parameter sampling distributions', r'$dp$ [hPa]', r'$r_p$ [km]', r'$|\phi_c|$ [°]', r'$q$ factor'
    ]
    ltx_df = ltx_df.reindex(new_order)

    # Display and Export
    print(ltx_df, '\n\n')
    
    # Use standard tabular with paragraph 'p' columns instead of tabularx 'X'
    latex_table = ltx_df.to_latex(
        escape=False,
        index=True,
        column_format='p{3cm} | p{2.2cm} p{2.2cm} p{2.2cm} p{2.2cm} p{2.2cm}',
        caption="Vortex perturbation parameters for each experiment.",
        label="tab:vortex_params"
    )

    # Note: We no longer string-replace \begin{tabular} with tabularx. 
    print(latex_table, '\n\n')

    final_output = (
        r"\renewcommand{\arraystretch}{1.5}" + "\n" # Lowered slightly as \newline creates real vertical space
        + latex_table
    )
    print(final_output)

if __name__ == '__main__':
    main()
