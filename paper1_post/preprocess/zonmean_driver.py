import sys
import os
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
from paths import ARCHRT, CAMGR
import subprocess
import signal

VARS = 'Z3,T,Q'
CONS = True

for ii, ar in enumerate(ARCHRT):
   hpth = os.path.join(ar, 'atm/hist/')

   proc = subprocess.Popen(f"qcmd -q casper -l walltime=02:21:00 -l select=1:ncpus=16:mem=256GB -A UCIS0005 python3 -u ux_zonmean.py\
                    {hpth} {CAMGR} {str(CONS)} {VARS} &> zmdriver.out{ii}", shell=True)
