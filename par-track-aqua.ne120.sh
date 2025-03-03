#!/bin/bash -l

##==============================================================
#PBS -N tempest.par.250225
#PBS -A UPSU0063
#PBS -l walltime=01:00:00
#PBS -q main
#PBS -j oe
#PBS -k eod
#PBS -l select=3:ncpus=128:mpiprocs=128
#PBS -m abe
#PBS -M jvp5930@psu.edu
################################################################

TEMPESTEXTREMESDIR=/glade/work/zarzycki/derecho/tempestextremes/
###TEMPESTEXTREMESDIR="/glade/u/home/jpan/miniconda3/pkgs/tempest-extremes-2.2.2-nompi_h2cd73b2_100/"

UQSTR=b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250225_unseed_all
PATHTOFILES=/glade/derecho/scratch/jpan/archive/${UQSTR}/atm/hist/
CONNECTDAT=/glade/u/home/jpan/ne120np4_connect_v2.dat
###CONNECTDAT=/glade/u/home/zarzycki/tempest-scripts/for-ben/ne120.connect_v2.dat
CONNECTFLAG=""
TOPOFILE=""

############ TRACKER MECHANICS #####################

cd /glade/u/home/jpan/aquaptc/tempest

###ls $PATHTOFILES

DATESTRING=`date +"%s%N"`
FILELISTNAME=filelist.txt.${DATESTRING}
TRAJFILENAME=trajectories.txt.${UQSTR}
touch $FILELISTNAME

ignoreyear=999999999
FILES=$(ls "${PATHTOFILES}"/*.h1i.*.nc | grep -v "$ignoreyear-")
###echo $FILES
for f in $FILES
do
  #echo "${f};${TOPOFILE}" >> $FILELISTNAME
  echo "${f}" >> $FILELISTNAME
done

starttime=$(date -u +"%s")

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

STR_DETECT="--verbosity 0 --timestride 2 --in_connect ${CONNECTDAT} --out cyclones_tempest.${DATESTRING} --closedcontourcmd PSL,${DCU_PSLFOMAG},${DCU_PSLFODIST},0;_DIFF(Z300,Z500),${DCU_WCFOMAG},${DCU_WCFODIST},${DCU_WCMAXOFFSET} --mergedist ${DCU_MERGEDIST} --searchbymin PSL --outputcmd PSL,min,0;_VECMAG(UBOT,VBOT),max,2"
echo "calling mpiexec"
mpiexec $TEMPESTEXTREMESDIR/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STR_DETECT} </dev/null

cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
#rm cyclones_tempest.${DATESTRING}*

pwd
${TEMPESTEXTREMESDIR}/bin/StitchNodes --format "ncol,lon,lat,slp,wind" --allow_repeated_times --range ${SN_TRAJRANGE} --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --threshold "wind,>=,${SN_MINWIND},${SN_MINLEN};lat,<=,${SN_MAXLAT},${SN_MINLEN};lat,>=,-${SN_MAXLAT},${SN_MINLEN}"

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

printf "${tottime}\n" >> timing.txt

#rm -v cyclones.${DATESTRING}.dat
#rm -v ${FILELISTNAME}
#rm -v log*txt
