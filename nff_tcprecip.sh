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

SPTH=17
GCD=8
UQSTR=b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl
PATHTOFILES=/glade/derecho/scratch/jpan/archive/${UQSTR}/atm/hist_0010_h1i/
DIRO=/glade/derecho/scratch/jpan/archive/${UQSTR}/atm/nff_${SPTH}mps
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
###FILES=$(ls "${PATHTOFILES}"/*.h1i.*.nc | grep -v "$ignoreyear-")
FILES=$(find "${PATHTOFILES}" -maxdepth 1 -name "*.h1i.*.nc" ! -name "*$ignoreyear-*" | sort)
###echo $FILES
for f in $FILES
do
  #echo "${f};${TOPOFILE}" >> $FILELISTNAME
  echo "${f}" >> $FILELISTNAME
  FILENAME=$(basename "$f")
  echo "${DIRO}/${FILENAME}.nff_${SPTH}mps" >> $OUTLISTNAME
done

starttime=$(date -u +"%s")

BINW=0.25
BINC=$( echo "$GCD / $BINW + 1" | bc )

UMEANRMV="_DIFF(UBOT,_MEAN{${GCD}}(UBOT))"
VMEANRMV="_DIFF(VBOT,_MEAN{${GCD}}(VBOT))"

STR_NFE1="--in_nodefile ${TRAJFILENAME} --in_nodefile_type SN --in_fmt ${SN_FMT} --in_data_list ${FILELISTNAME} --in_connect ${CONNECTDAT} --out_nodefile ${TRAJFILENAME}.radspd --out_fmt ${SN_FMT},radspd,r${SPTH} --calculate radspd=radial_wind_profile(${UMEANRMV},${VMEANRMV},${BINC},${BINW});r${SPTH}=lastwhere(radspd,>=,${SPTH})*${BINW}"

$TEMPESTEXTREMESDIR/bin/NodeFileEditor ${STR_NFE1}

###STR_NFE2="--in_nodefile ${TRAJFILENAME}.radspd --in_nodefile_type SN --in_fmt ${SN_FMT},radspd --in_data_list ${FILELISTNAME} --in_connect ${CONNECTDAT} --out_nodefile ${TRAJFILENAME}.radspd --out_fmt ${SN_FMT},radspd,r8 --calculate r8=lastwhere(radspd,>=,8)*${BINW}"

###$TEMPESTEXTREMESDIR/bin/NodeFileEditor ${STR_NFE2}

STR_NFF="--in_nodefile ${TRAJFILENAME}.radspd --in_nodefile_type SN --in_fmt ${SN_FMT},radspd,r${SPTH} --in_data_list ${FILELISTNAME} --in_connect ${CONNECTDAT} --out_data_list ${OUTLISTNAME} --maskvar TC_R${SPTH} --bydist r${SPTH}"
###--bycontour _PROD(_SIGN(lat),_CURL{8,1.0}(U850,V850)),-1e-5,5.5,0.5"
###PSL,${DCU_PSLFOMAG},${DCU_PSLFODIST},0"

###echo "calling mpiexec"
###mpiexec --display-allocation --display-map --report-bindings -n 16 $TEMPESTEXTREMESDIR/bin/NodeFileFilter ${STR_NFF} </dev/null
$TEMPESTEXTREMESDIR/bin/NodeFileFilter ${STR_NFF}
