import sys
import os
sys.path.append('/glade/u/home/jpan/aquaptc/bmom6_tcdiag/paper1_post/')
from paths import ARCHRT
import subprocess

filtvars = 'TAUX:PRECT:V850:Q850' #need to regrid TC masks to ocean to filter evap
rgcd = 8
wspd = 4
invert = False

for ar in ARCHRT:
   casenm = ar.rstrip('/').split('/')[-1]
   atmpth = os.path.join(ar, 'atm/hist')
   #subprocess.run(f"qsub -v 'UQSTR={casenm},PATHTOFILES={atmpth},FILTVARS={filtvars},GCD={rgcd},SPTH={wspd}' nff_gen_TC_masks.sh", shell=True, check=True)
   subprocess.Popen(f"UQSTR='{casenm}' PATHTOFILES='{atmpth}' FILTVARS='{filtvars}' GCD={rgcd} SPTH={wspd} INVERT={str(invert).lower()} ./nff_gen_TC_masks.sh &> nffout.{casenm}{'_invert' if invert else ''}", shell=True)
