import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FILI = ['250415_unseed.parquet', '250417_ctrl.parquet', '250416_seed1x1.parquet']
ALIASES = ['UNSEED', 'CTRL', 'SEED']
PBINS = np.arange(830, 1030, 5)

dfs = [pd.read_parquet(f) for f in FILI]
pmins = [df.groupby('stmnum')['pres'].min() / 100 for df in dfs]
meanvals = [pm.mean() for pm in pmins]
print(min(pmins[2]))

plt.rc('font', size=17)
plt.rcParams['figure.figsize'] = (24, 5)
fig, axes = plt.subplots(1, 3, sharex=True, sharey=True)
[axes[ii].hist(pm, edgecolor='black', bins=PBINS) for ii, pm in enumerate(pmins)]
axes[0].set_ylabel('# of storms per bin')
[ax.set_xlabel('Minimum SLP [hPa]') for ax in axes]
[ax.set_title(ALIASES[ii] + ' (mean = %.1f)' % meanvals[ii]) for ii, ax in enumerate(axes)]
[ax.axvline(meanvals[ii], linewidth=2.5, color='black', linestyle='dashed') for ii, ax in enumerate(axes)]
[ax.tick_params(axis='x', labelrotation=45) for ax in axes]
[axes[ii].set_title(ss, loc='left') for ii, ss in enumerate(['(a)', '(b)', '(c)'])]
axes[0].set_xticks(np.arange(840, 1030, 20))
plt.savefig('lmi_psl.svg', bbox_inches='tight')
plt.show()
