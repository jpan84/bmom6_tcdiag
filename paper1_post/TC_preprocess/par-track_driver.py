import sys
import os
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
from paths import ARCHRT
import subprocess

for ar in ARCHRT:
   casenm = ar.rstrip('/').split('/')[-1]
   atmpth = os.path.join(ar, 'atm/hist')
   subprocess.run(f"qsub -v UQSTR='{casenm}' -v PATHTOFILES='{atmpth}' par-track-aqua.derecho_worker.sh", shell=True, check=True)
