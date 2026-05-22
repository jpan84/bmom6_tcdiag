import sys
import os
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
from paths import ARCHRT
import subprocess

#EXPDIR = '/glade/derecho/scratch/jpan/archive/%s/atm'
#EXPS = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250702_unseed2hPa6m', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']

for ar in ARCHRT:
   atmpth = os.path.join(ar, 'atm/')
   subprocess.run(f"qsub -v EXP_DIR='{atmpth}' hy2pres_gnupar.sh", shell=True, check=True)
