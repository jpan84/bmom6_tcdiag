#!/bin/bash

CONNECTDAT="/glade/u/home/jpan/ne120np4_connect_v2.dat"
TEMPESTFILE="./unseed_dn_tempest"
TEMPESTEXTREMESDIR="/glade/work/zarzycki/derecho/tempestextremes_noMPI"
TESTFILE="/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl2e-5/rest/0005-02-14-00000/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl2e-5.cam.h1i.0005-02-14-00000.nc"
TESTFILE="/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl1deg/rest/0005-02-14-00000/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250303_unseed_curl1deg.cam.h1i.0005-02-14-00000.nc"

# Do TE to find candidates in final timeslice
###STR_DETECT="--verbosity 3 --in_connect ${CONNECTDAT} --out ${TEMPESTFILE} --closedcontourcmd PSL,400.0,5.5,0;_DIFF(Z300,Z500),-15,6.5,0.5;_PROD(_SIGN(lat),_CURL{8,2.0}(U850,V850)),-2e-5,5.5,0.5;T500,-2.0,6.5,0.5 --minlat -30.0 --maxlat 30.0 --mergedist 6.0 --searchbymin PSL --outputcmd PSL,min,0"
###${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data ${TESTFILE} ${STR_DETECT} >> "${TEMPESTFILE}.log"

# Do TE to find candidates in final timeslice
STR_DETECT="--verbosity 3 --in_connect ${CONNECTDAT} --out ${TEMPESTFILE} --closedcontourcmd PSL,400.0,5.5,0;_DIFF(Z300,Z500),-15,6.5,0.5 --minlat -30.0 --maxlat 30.0 --mergedist 6.0 --searchbymin PSL --outputcmd PSL,min,0;_PROD(_SIGN(lat),_CURL{8,1.0}(U850,V850)),max,0.5"
${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data ${TESTFILE} ${STR_DETECT} > "${TEMPESTFILE}.log"
