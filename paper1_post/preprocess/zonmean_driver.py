import sys
import os
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
from paths import ARCHRT, CAMGR
import subprocess
import signal

VARS = 'Z3,T,Q'
CONS = True

for ii, ar in enumerate(ARCHRT):
   atmpth = os.path.join(ar, 'atm/')

   proc = subprocess.Popen(f"setsid nohup qcmd -q casper -l walltime=02:21:00 -l select=1:ncpus=16:mem=256GB -A UCIS0005 python3 -u ux_zonmean.py\
                    {atmpth} {CAMGR} {str(CONS)} {VARS} &> zmdriver.out{ii}", shell=True)
