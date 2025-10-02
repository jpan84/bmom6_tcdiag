import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FILI = ['250415_unseed.parquet', '250417_ctrl.parquet', '250416_seed1x1.parquet']
ALIASES = ['UNSEED', 'CTRL', 'SEED']
PBINS = np.arange(830, 1030, 10)

dfs = [pd.read_parquet(f) for f in FILI]
pmins = [df.groupby('stmnum')['pres'].min() / 100 for df in dfs]
print(min(pmins[2]))

plt.rc('font', size=16)
plt.rcParams['figure.figsize'] = (30, 6)
fig, axes = plt.subplots(1, 3, sharex=True, sharey=True)
[axes[ii].hist(pm, edgecolor='black', bins=PBINS) for ii, pm in enumerate(pmins)]
axes[0].set_ylabel('# of storms per bin')
[ax.set_xlabel('Minimum SLP [hPa]') for ax in axes]
[ax.set_title(ALIASES[ii]) for ii, ax in enumerate(axes)]
[axes[ii].set_title(ss, loc='left') for ii, ss in enumerate(['(a)', '(b)', '(c)'])]
axes[0].set_xticks(np.arange(840, 1030, 20))
plt.savefig('lmi_psl.png', bbox_inches='tight')
plt.show()
