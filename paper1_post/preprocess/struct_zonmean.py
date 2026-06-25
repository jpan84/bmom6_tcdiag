#cdo zonmean -cat -apply,-selvar,thetao [ b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250702_uns
#eed2hPa6m.mom6.hmcust_avg.0009*.nc ] testcdozm.nc

import os
import sys
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
from paths import ARCHRT, OCN_CUST
import subprocess

#DIRI = sys.argv[1]
OHST = 'ocn/hist/'
TAPE0, TAPE1 = '*hmcust_avg*.nc', '*hm.00*.nc' #file names to process
VARS = sys.argv[1] #example of comma-separated (no space): PS,PRECT,QFLX

#subprocess.run('module load cdo', shell=True, check=True)
for ii, ar in enumerate(ARCHRT):
   hfil = os.path.join(ar, OHST, TAPE0 if OCN_CUST[ii] else TAPE1)
   filo = os.path.join(ar, OHST, 'cdo_zm_' + '_'.join(VARS.split(',')) + '.nc')

   proc = None
   if OCN_CUST[ii]:
      proc = subprocess.Popen(f"module load cdo && qcmd -q casper -l walltime=02:20:00 -l select=1:ncpus=8:mem=128GB -A UCIS0005 -- \
                                cdo zonmean -cat '-apply,-selvar,{VARS}' [ {hfil} ] {filo} &> cdozm.out{ii}", shell=True)
   else:
      proc = subprocess.Popen(f"module load cdo && qcmd -q casper -l walltime=00:40:00 -l select=1:ncpus=8:mem=128GB -A UCIS0005 -- \
                                cdo zonmean -selvar,{VARS} -cat [ {hfil} ] {filo} &> cdozm.out{ii}", shell=True)
   print(proc.args)
