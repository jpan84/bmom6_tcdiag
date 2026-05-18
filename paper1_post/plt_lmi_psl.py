import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DIRI = '~/aquaptc/tempest/'
FILI = ['250702_unseed2hPa6m.parquet', '250415_unseed.parquet', '250417_ctrl.parquet', '251229_seedmatch.parquet', '250416_seed1x1.parquet']
ALIASES = ['UNSEED_90', 'UNSEED_50', 'CTRL', 'SEED_50', 'SEED_150']
FILI = ['250415_unseed.parquet', '250417_ctrl.parquet', '251229_seedmatch.parquet','250702_unseed2hPa6m.parquet', '250416_seed1x1.parquet']
ALIASES = ['UNSEED_50', 'CTRL', 'SEED_50', 'UNSEED_90', 'SEED_150']
PBINS = np.arange(830, 1030, 5)

SZNS = ['DJF', 'MAM', 'JJA', 'SON']
SZMOs = [{12, 1, 2}, {3, 4, 5}, {6, 7, 8}, {9, 10, 11}]
WARM = [-1, -1, 1, 1] #sign of warm hemi

dfs = [pd.read_parquet(os.path.join(DIRI, f)) for f in FILI]
totyrs = [(df.index[-1] - df.index[0]).days / 365. for df in dfs]
dfs_filtered = [
    df[np.logical_or.reduce([
        (df.index.month.isin(mos) & (np.sign(df['lat']) == sgn)) 
        for mos, sgn in zip(SZMOs, WARM)
    ])] 
    for df in dfs
]
pmins = [df.groupby('stmnum')['pres'].min() / 100 for df in dfs_filtered]
meanvals = [pm.mean() for pm in pmins]
print(min(pmins[2]))

plt.rc('font', size=17)
plt.rcParams['figure.figsize'] = (24, 10)
fig, axes = plt.subplots(2, 3, sharex=True, sharey=True)
useaxes = [ax for ax in axes.ravel()]
useaxes[-2].set_axis_off()
useaxes.pop(-2)
[useaxes[ii].hist(pm, edgecolor='black', bins=PBINS, weights=np.ones_like(pm) / totyrs[ii]) for ii, pm in enumerate(pmins)]
useaxes[0].set_ylabel('# of warm-season TCs per year')
[ax.set_xlabel('Minimum SLP [hPa]') for ax in useaxes]
[ax.tick_params(labelbottom=True, labelleft=True, top=True, right=True) for ax in useaxes]
[ax.set_title(ALIASES[ii] + ' (mean = %.1f)' % meanvals[ii]) for ii, ax in enumerate(useaxes)]
[ax.axvline(meanvals[ii], linewidth=2.5, color='black', linestyle='dashed') for ii, ax in enumerate(useaxes)]
[ax.tick_params(axis='x', labelrotation=45) for ax in useaxes]
[useaxes[ii].set_title(ss, loc='left') for ii, ss in enumerate(['(a)', '(b)', '(c)', '(d)', '(e)'])]
useaxes[0].set_xticks(np.arange(840, 1030, 20))
fig.tight_layout()
plt.savefig('lmi_psl.svg', bbox_inches='tight')
plt.show()
