import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PQ1 = sys.argv[1]
PQ2 = sys.argv[2]
VARS = ['clat', 'psmin', 'rp', 'dp', 'zp', 'exppr']
bn = lambda x: os.path.basename(x)

df1 = pd.read_parquet(PQ1).dropna()
df2 = pd.read_parquet(PQ2).dropna()

for vv in VARS:
   cat = np.concatenate((df1[vv], df2[vv]))
   ct, be = np.histogram(cat, bins='auto')

   plt.hist(df1[vv], bins=be, alpha=1, density=True)
   plt.hist(df2[vv], bins=be, alpha=0.3, density=True)
   plt.title(vv)

   plt.savefig(os.path.join(os.path.dirname(PQ1), bn(PQ1).split('_')[0] + '_' + bn(PQ2).split('_')[0] + '_' + str(vv) + '.png'))
   plt.close()
