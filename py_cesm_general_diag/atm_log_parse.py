#construct time series of CAM energy/mass checks from concatenated archived logs
import gzip
import os
#import glob
#import re
import numpy as np
import matplotlib.pyplot as plt

ALIS = ['250415_unseed', '250417_ctrl', '250416_seed1x1']
CLR = ['blue', 'orange', 'green']
FILI = '/glade/u/home/jpan/work/MOM6_CASEDIRS/250708_logs_backup/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/atm.log.cat.gz'

def main():
   for ii, al in enumerate(ALIS):
      lines = []
      with gzip.open(FILI % al, 'rt', encoding='utf-8') as f:
         lines = f.readlines()
   
      spl = [ln.split() for ln in lines]
      te_lns = []
      for ln in spl:
         if len(ln) >= 2 and ln[0] == 'nstep,' and ln[1] == 'te':
            te_lns.append([float(st) for st in ln[2:]])
   
      te_arr = np.array(te_lns)
      #print(te_arr)
      te_arr[:, 1:3] *= 4 * np.pi * (6.371e6)**2
   
      plt.plot(te_arr[:, 0] / 192 / 365 + 1, te_arr[:, 1], color=CLR[ii], label=al, linewidth=0.8)
      #plt.plot(te_arr[:, 0] / 192 / 365 + 1, float(te_arr[:, 1].mean()), color=CLR[ii], linestyle='--', linewidth=0.5)
      plt.hlines(te_arr[:, 1].mean(), 5, 16, color=CLR[ii], linestyle='--', linewidth=0.5)

   plt.rc('font', size=16)
   plt.xlabel('year')
   plt.ylabel('CAM atm.log total energy [J]')
   plt.legend()
   plt.show()

if __name__ == '__main__':
   main()
