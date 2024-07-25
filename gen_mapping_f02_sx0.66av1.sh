##qcmd -A UPSU0063 -- /glade/u/home/jpan/work/CESM2.3_alpha17b/cime/tools/mapping/gen_mapping_files/gen_cesm_maps.sh \
##     --fileocn  /glade/work/youweima/f09_sx0.66av1/sx0.66av1_grid_SCRIP.nc\
##     --fileatm  /glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.23x0.31_141008.nc  \
##     --nameocn  sx0.66av1 \
##     --nameatm  fv0.23x0.31 \
##     --serial

##qsub -- ~/VR-CESM-toolkit/gen_mapping/merged-mapping.sh \
##    --atmName fv0.23x0.31 \
##    --ocnName sx0.66av1 \
##    --atmGridName /glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.23x0.31_141008.nc \
##    --ocnGridName /glade/work/youweima/f09_sx0.66av1/sx0.66av1_grid_SCRIP.nc \
##    --wgtFileDir /glade/work/jpan/f02_sx0.66av1

module load ncl esmf/8.6.0

nclscr="/glade/u/home/jpan/VR-CESM-toolkit/gen_mapping/gen_X_to_Y_wgts.ncl"
atmName="fv0.23x0.31"
atmGridName="/glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.23x0.31_141008.nc"
ocnName="sx0.66av1"
ocnGridName="/glade/work/jpan/f09_sx0.66av1/sx0.66av1_grid_SCRIP.nc"
ocnMeshName="/glade/work/jpan/f09_sx0.66av1/ESMF_mesh_sx0.66av1.nc"
wgtFileDir="/glade/work/jpan/f02_sx0.66av1"
cdate=`date +%y%m%d`

###interp_method="conserve"
###ncl ${nclscr} 'srcName="'${atmName}'"' 'srcGridName="'${atmGridName}'"' 'dstName="'${ocnName}'"' 'dstGridName="'${ocnGridName}'"' 'wgtFileDir="'${wgtFileDir}'"' 'InterpMethod="'${interp_method}'"'

# do ATM2OCN_SMAPNAME and ATM2OCN_VMAPNAME (blin)
interp_method="bilinear"   # bilinear, patch, conserve
ncl ${nclscr} 'srcName="'${atmName}'"' 'srcGridName="'${atmGridName}'"' 'dstName="'${ocnName}'"' 'dstGridName="'${ocnGridName}'"' 'wgtFileDir="'${wgtFileDir}'"' 'InterpMethod="'${interp_method}'"'

# do OCN2ATM_FMAPNAME and OCN2ATM_SMAPNAME (aave)
###interp_method="conserve"   # bilinear, patch, conserve
###ncl ${nclscr} 'srcName="'${ocnName}'"' 'srcGridName="'${ocnGridName}'"' 'dstName="'${atmName}'"' 'dstGridName="'${atmGridName}'"' 'wgtFileDir="'${wgtFileDir}'"' 'InterpMethod="'${interp_method}'"'

###ncremap -s ${atmGridName} -g ${ocnGridName} -m ${wgtFileDir}/map_${atmName}_TO_${ocnName}_aave.nc -a aave

ESMF_RegridWeightGen -s ${atmGridName} -d ${ocnMeshName} -w ${wgtFileDir}/map_${atmName}_TO_${ocnName}_aave_${cdate}.nc --method conserve --ignore_unmapped --ignore_degenerate
ESMF_RegridWeightGen -s ${ocnMeshName} -d ${atmGridName} -w ${wgtFileDir}/map_${ocnName}_TO_${atmName}_aave_${cdate}.nc --method conserve --ignore_unmapped --ignore_degenerate
