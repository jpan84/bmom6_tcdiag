#!/bin/bash
####PBS -N tempest.NFF
####PBS -A UCIS0005
####PBS -l select=1:ncpus=18:mpiprocs=18:mem=80GB
####PBS -l walltime=0:05:00
####PBS -q casper@casper-pbs
####PBS -j oe
####PBS -m a
####PBS -M jvp5930@psu.edu
###################################################################
###
###starttime=$(date -u +"%s")
###
#### Casper modules
###module load intel
###module load openmpi
###module load parallel
###module load nco
###
###echo "MPI info"
###nproc
###cat /proc/cpuinfo | grep processor | wc -l
###grep "cpu cores" /proc/cpuinfo | uniq
###grep -c processor /proc/cpuinfo

###TEMPESTEXTREMESDIR=/glade/work/zarzycki/derecho/tempestextremes/
TEMPESTEXTREMESDIR=/glade/work/zarzycki/tempestextremes_noMPI


UQSTR=b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed
PATHTOFILES=/glade/derecho/scratch/jpan/archive/${UQSTR}/atm/hist_0010_h1i/
DIRO=/glade/derecho/scratch/jpan/archive/${UQSTR}/atm/nff_tcprec
mkdir -p $DIRO
CONNECTDAT=/glade/u/home/jpan/ne120np4_connect_v2.dat
CONNECTFLAG=""
TOPOFILE=""

SN_FMT="lon,lat,slp,wind"

############ TRACKER MECHANICS #####################

###cd /glade/u/home/zarzycki/tempest-scripts/for-joshua
cd /glade/u/home/jpan/aquaptc/tempest

###ls $PATHTOFILES

DATESTRING=`date +"%s%N"`
FILELISTNAME=filelist.txt.${DATESTRING}
OUTLISTNAME=outlist.txt.${DATESTRING}
TRAJFILENAME=trajectories.txt.${UQSTR}
touch $FILELISTNAME $OUTLISTNAME

ignoreyear=999999999
FILES=$(ls "${PATHTOFILES}"/*.h1i.*.nc | grep -v "$ignoreyear-")
###echo $FILES
for f in $FILES
do
  #echo "${f};${TOPOFILE}" >> $FILELISTNAME
  echo "${f}" >> $FILELISTNAME
  FILENAME=$(basename "$f")
  echo "${DIRO}/${FILENAME}" >> $OUTLISTNAME
done

starttime=$(date -u +"%s")

DCU_PSLFOMAG=200.0
DCU_PSLFODIST=5.5

STR_NFF="--in_nodefile ${TRAJFILENAME} --in_nodefile_type SN --in_fmt ${SN_FMT} --in_data_list ${FILELISTNAME} --in_connect ${CONNECTDAT} --out_data_list ${OUTLISTNAME} --var PRECT,OMEGA500,OMEGA850,V850,T850,V500,T500,TAUX --maskvar TC8GCD5MPS --nearbyblobs _VECMAG(UBOT,VBOT),0.0,>,5.0,8.0"
###--bycontour _PROD(_SIGN(lat),_CURL{8,1.0}(U850,V850)),-1e-5,5.5,0.5"
###PSL,${DCU_PSLFOMAG},${DCU_PSLFODIST},0"

###echo "calling mpiexec"
###mpiexec --display-allocation --display-map --report-bindings -n 16 $TEMPESTEXTREMESDIR/bin/NodeFileFilter ${STR_NFF} </dev/null
$TEMPESTEXTREMESDIR/bin/NodeFileFilter ${STR_NFF}











exit
##############################################################

###orig
DCU_PSLFOMAG=200.0
DCU_PSLFODIST=5.5
DCU_WCFOMAG=-6.0
DCU_WCFODIST=6.5
DCU_WCMAXOFFSET=1.0
DCU_MERGEDIST=6.0
SN_TRAJRANGE=6.0
SN_TRAJMINLENGTH=10
SN_TRAJMAXGAP=3
SN_MAXTOPO=150.0
SN_MAXLAT=50.0
SN_MINWIND=10.0
SN_MINLEN=10

touch cyclones.${DATESTRING}

STR_DETECT="--verbosity 0 --timestride 6hr --in_connect ${CONNECTDAT} --out cyclones_tempest.${DATESTRING} --closedcontourcmd PSL,${DCU_PSLFOMAG},${DCU_PSLFODIST},0;_DIFF(Z300,Z500),${DCU_WCFOMAG},${DCU_WCFODIST},${DCU_WCMAXOFFSET};_PROD(_SIGN(lat),_CURL{8,1.0}(U850,V850)),-8e-5,5.5,0.5 --mergedist ${DCU_MERGEDIST} --searchbymin PSL --outputcmd PSL,min,0;_VECMAG(UBOT,VBOT),max,2"
echo "calling mpiexec"
mpiexec --display-allocation --display-map --report-bindings -n 16 $TEMPESTEXTREMESDIR/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STR_DETECT} </dev/null

cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
#rm cyclones_tempest.${DATESTRING}*

pwd
${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_connect ${CONNECTDAT} --in_fmt ${DN_FMT} --allow_repeated_times --range ${SN_TRAJRANGE} --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --threshold "wind,>=,${SN_MINWIND},${SN_MINLEN};lat,<=,${SN_MAXLAT},${SN_MINLEN};lat,>=,-${SN_MAXLAT},${SN_MINLEN}"

NFEMEANRAD=8.0

###${TEMPESTEXTREMESDIR}/bin/NodeFileEditor --in_nodefile ${TRAJFILENAME} --in_nodefile_type "SN" --in_data_list "${FILELISTNAME}" --in_connect ${CONNECTDAT} --in_fmt "${DN_FMT}" --out_fmt ${NFE_FMT} --calculate "prect=_MEAN{${NFEMEANRAD}}(PRECT);olr=_MEAN{${NFEMEANRAD}}(FLUT);sst=_MEAN{${NFEMEANRAD}}(SST);ace0=eval_ace(UBOT,VBOT,0.25);radialpsl=radial_profile(PSL,32,0.25)" --out_nodefile "${TRAJFILENAME}.NFE"

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

printf "${tottime}\n" >> timing.txt

rm -v cyclones.${DATESTRING}.dat
rm -v ${FILELISTNAME}
rm -v log*txt
