import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FILI = ['250415_unseed.kde.pickle', '250417_ctrl.kde.pickle', '250416_seed1x1.kde.pickle']
PQTS = ['250415_unseed.parquet', '250417_ctrl.parquet', '250416_seed1x1.parquet']

pkls = []
for fl in FILI:
   with open(fl, 'rb') as ofl:
      pkls.append(pickle.load(ofl))
#print(pkls)

xgr, ygr = np.mgrid[1:365:365j, -25:25:100j]
print(xgr[:, 0], ygr[0, ygr.shape[1]//2:], ygr[0, :ygr.shape[1]//2])

kdeobjs = [pk[('genday', 'genlat')][0] for pk in pkls]
nstms = [pk[('genday', 'genlat')][1] for pk in pkls]
print(kdeobjs[0])
print(nstms)

zvals = [ko.evaluate(np.vstack([xgr.ravel(), ygr.ravel()])).reshape(xgr.shape) * nstms[ii] for ii, ko in enumerate(kdeobjs)]
print(zvals[0])
znh = [zv[:, ygr.shape[1] // 2:] for zv in zvals]
zsh = [np.roll(zv[:, :ygr.shape[1] // 2][:, ::-1], 182, axis=0) for zv in zvals]
zplt = [znh[ii] + zsh[ii] for ii in range(len(znh))]

plt.contourf(xgr[:, ygr.shape[1] // 2:], ygr[:, ygr.shape[1] // 2:], zsh[1])
plt.show()

dfs = [pd.read_parquet(pq) for pq in PQTS]
print(dfs[0])
