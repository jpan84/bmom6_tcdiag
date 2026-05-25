import sys
import os
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
from paths import ARCHRT
import subprocess

for ar in ARCHRT:
   atmpth = os.path.join(ar, 'atm/')
   subprocess.run(f"qsub -v EXP_DIR='{atmpth}' hy2pres_gnupar.sh", shell=True, check=True)
