#!/bin/bash

TEMPESTEXTREMESDIR=/glade/work/zarzycki/tempestextremes_noMPI

###TODO: pass args from driver (SPTH, GCD, UQSTR, PATHTOFILES, FILTVARS)
###SPTH=4
###GCD=8
###UQSTR=b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl
###PATHTOFILES=/glade/derecho/scratch/jpan/jpan_tcfields/${UQSTR}/hist_0012-0014_h1i/
DIRO=${PATHTOFILES}/../nff_${SPTH}mps
mkdir -p $DIRO
CONNECTDAT=/glade/u/home/jpan/ne120np4_connect_v2.dat
CONNECTFLAG=""
TOPOFILE=""


SN_FMT="lon,lat,slp,wind"

############ TRACKER MECHANICS #####################

DATESTRING=`date +"%s%N"`
FILELISTNAME=filelist.txt.${DATESTRING}
OUTLISTNAME=outlist.txt.${DATESTRING}
TRAJFILENAME=trajectories.txt.${UQSTR}
touch $FILELISTNAME $OUTLISTNAME

ignoreyear=999999999
DATE_LIM="0007-02"
###FILES=$(ls "${PATHTOFILES}"/*.h1i.*.nc | grep -v "$ignoreyear-")
FILES=$(find "${PATHTOFILES}" -maxdepth 1 -name "*.h1i.*.nc" ! -name "*$ignoreyear-*" | sort)
###echo $FILES
for f in $FILES
do
  FILE_DATE=$(echo "$f" | grep -oE '[0-9]{4}-[0-9]{2}')

  if [[ "$FILE_DATE" < "${DATE_LIM}" ]]; then
    echo "${f}" >> $FILELISTNAME
    FILENAME=$(basename "$f")
    echo "${DIRO}/${FILENAME}.nff_${SPTH}mps" >> $OUTLISTNAME
  fi

done

starttime=$(date -u +"%s")

BINW=0.25
BINC=$( echo "$GCD / $BINW + 1" | bc )

STR_NFE1="--in_nodefile ${TRAJFILENAME} --in_nodefile_type SN --in_fmt ${SN_FMT} --in_data_list ${FILELISTNAME} --in_connect ${CONNECTDAT} --out_nodefile ${TRAJFILENAME}.radspd${SPTH} --out_fmt ${SN_FMT},radspd,r${SPTH} --calculate radspd=radial_wind_profile(UBOT,VBOT,${BINC},${BINW});r${SPTH}=lastwhere(radspd,>=,${SPTH})*${BINW}"

$TEMPESTEXTREMESDIR/bin/NodeFileEditor ${STR_NFE1}

STR_NFF="--in_nodefile ${TRAJFILENAME}.radspd${SPTH} --in_nodefile_type SN --in_fmt ${SN_FMT},radspd,r${SPTH} --in_data_list ${FILELISTNAME} --in_connect ${CONNECTDAT} --out_data_list ${OUTLISTNAME} --maskvar TC_R${SPTH} --var ${FILTVARS//:/,} --bydist r${SPTH}"

$TEMPESTEXTREMESDIR/bin/NodeFileFilter ${STR_NFF}
