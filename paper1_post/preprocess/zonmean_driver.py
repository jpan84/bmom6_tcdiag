import sys
import os
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
from paths import ARCHRT, CAMGR
import subprocess
import signal

#VARS = 'Z3,T,Q'
VARS = sys.argv[1]
TAPE = sys.argv[2] if len(sys.argv) > 2 else 'atm/hist/*h0a.00[0-9]*.nc' #path and file names to process
CONS = False

for ii, ar in enumerate(ARCHRT):
   hpth = os.path.join(ar, TAPE)
   #esc_hpth = hpth.replace('*', r'\*').replace('[', r'\[').replace(']', r'\]') #attempted to prevent bash from expanding wildcards
   pthptr = f'zmdriver.pthptr{ii}'
   with open(pthptr, 'w') as f:
      f.write(hpth)

   proc = subprocess.Popen(f"qcmd -q casper -l walltime=03:00:00 -l select=1:ncpus=16:mem=256GB -A UCIS0005 python3 -u ux_zonmean.py\
                    {pthptr} {CAMGR} {str(CONS)} {VARS} &> zmdriver.out{ii}", shell=True)
   print(proc.args)
