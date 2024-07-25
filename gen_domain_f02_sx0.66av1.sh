cd $WORK/f02_sx0.66av1

qcmd -A UPSU0063 -- /glade/u/home/youweima/work/CESM2.3/cime/tools/mapping/gen_domain_files/gen_domain \
           -m ./map_sx0.66av1_TO_fv0.23x0.31_aave_240725.nc \
           -o sx0.66av1 \
           -l fv0.23x0.31
