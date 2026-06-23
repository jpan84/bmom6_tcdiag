#!/bin/bash -l
#PBS -N tempest.CZ.worker
#PBS -A UPSU0032
#PBS -l select=1:ncpus=64:mpiprocs=64:mem=128GB
#PBS -l walltime=06:00:00
#PBS -q develop
#PBS -j oe
#PBS -m a
#PBS -M jvp5930@psu.edu
################################################################

starttime=$(date -u +"%s")
execdir=$(pwd)
echo "$starttime $PBS_NP"

set -e
# Derecho modules
module load cmake
module load ncarenv
module load ncarcompilers
module load intel
module load cray-mpich
module load netcdf
export LD_LIBRARY_PATH=${NCAR_LDFLAGS_NETCDF}64:${NCAR_LDFLAGS_NETCDF}:${LD_LIBRARY_PATH}

echo "MPI info"
nproc
cat /proc/cpuinfo | grep processor | wc -l
grep "cpu cores" /proc/cpuinfo | uniq
grep -c processor /proc/cpuinfo

TEMPESTEXTREMESDIR=/glade/work/zarzycki/derecho/tempestextremes/

###UQSTR="$1"
###b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch
###PATHTOFILES="$2"
###/glade/derecho/scratch/jpan/archive/${UQSTR}/atm/hist/
CONNECTDAT=/glade/u/home/jpan/ne120np4_connect_v2.dat
CONNECTFLAG=""
TOPOFILE=""

DN_FMT="lon,lat,slp,wind"

############ TRACKER MECHANICS #####################

DATESTRING=`date +"%s%N"`
TMPDIR=/glade/work/jpan/tmpTE.${DATESTRING}/
mkdir -p ${TMPDIR}
cd ${TMPDIR} || exit 1

cleanup() {
  rm -rf "${TMPDIR}"
}
trap cleanup EXIT

FILELISTNAME=filelist.txt.${DATESTRING}
TRAJFILENAME=${execdir}/trajectories.txt.${UQSTR}
touch $FILELISTNAME

DATE_LIM="0020-02"
for f in "${PATHTOFILES}"/*.h1i.*.nc; do
  # Extract the YYYY-MM part
  FILE_DATE=$(echo "$f" | grep -oE '[0-9]{4}-[0-9]{2}')

  # If the date is strictly older than DATE_LIM, process it.
  if [[ "$FILE_DATE" < "${DATE_LIM}" ]]; then
    echo "${f}" >> "$FILELISTNAME"
  fi
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

STR_DETECT="--verbosity 0 --timestride 6hr --in_connect ${CONNECTDAT} --out cyclones_tempest.${DATESTRING} --closedcontourcmd PSL,${DCU_PSLFOMAG},${DCU_PSLFODIST},0;_DIFF(Z300,Z500),${DCU_WCFOMAG},${DCU_WCFODIST},${DCU_WCMAXOFFSET};_PROD(_SIGN(lat),_CURL{8,1.0}(U850,V850)),-8e-5,5.5,0.5 --mergedist ${DCU_MERGEDIST} --searchbymin PSL --outputcmd PSL,min,0;_VECMAG(UBOT,VBOT),max,2"
echo "calling mpiexec"
mpiexec  -n 64 $TEMPESTEXTREMESDIR/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STR_DETECT} </dev/null
###--display-allocation --display-map --report-bindings

cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
#rm cyclones_tempest.${DATESTRING}*

pwd
${TEMPESTEXTREMESDIR}/bin/StitchNodes --in_connect ${CONNECTDAT} --in_fmt ${DN_FMT} --allow_repeated_times --range ${SN_TRAJRANGE} --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --threshold "wind,>=,${SN_MINWIND},${SN_MINLEN};lat,<=,${SN_MAXLAT},${SN_MINLEN};lat,>=,-${SN_MAXLAT},${SN_MINLEN}"

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

printf "${tottime}\n" >> timing.txt

rm -rf "${TMPDIR}"
###rm -v cyclones.${DATESTRING}.dat
###rm -v ${FILELISTNAME}
###rm -v log*txt
