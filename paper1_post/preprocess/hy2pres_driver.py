import subprocess

EXPDIR = '/glade/derecho/scratch/jpan/archive/%s/atm'

EXPS = ['b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250702_unseed2hPa6m', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch', 'b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1']

for ex in EXPS:
   subprocess.run(f"qsub -v EXP_DIR='{EXPDIR % ex}' hy2pres_gnupar.sh", shell=True, check=True)
