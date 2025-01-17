#!/bin/bash

module load nco

CASENAME="b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.250113_seed3x3NH"
tapes=("rho2" "hm" "h")

for tape in "${tapes[@]}"; do
    echo $tape
    MOSTRECENT=`ls --color=never ./$CASENAME.mom6.$tape.*.nc -r | head -n 1`
    echo $MOSTRECENT
    LASTDATE=$(echo "$MOSTRECENT" | awk -F. '{print $(NF-1)}')
    echo $LASTDATE
    ###LASTYYYY=$(echo "$LASTDATE" | awk -F- '{print $(NF-2)}')
    ###echo $LASTYYYY
    ###LASTMM=$(echo "$LASTDATE" | awk -F- '{print $(NF-1)}')
    ###echo $LASTMM
    ###LASTDD=$(echo "$LASTDATE" | awk -F- '{print $NF}')
    ###echo $LASTDD
    IFS='-' read -r LASTYYYY LASTMM LASTDD <<< "$LASTDATE"
    echo $LASTDD

    ###noleap
    DDINMO=$(cal $LASTMM 0001 | awk 'NF {DAYS = $NF}; END {print DAYS}')
    echo $DDINMO

    if [[ "$LASTDD" == "$DDINMO" ]]; then
        MOFILES="$CASENAME.mom6.$tape.$LASTYYYY-$LASTMM-*.nc"
        ncra $MOFILES "$CASENAME.mom6.$tape.$LASTYYYY-$LASTMM.avg.nc"
        ###rm -v $MOFILES
        ###TODO: mv stuff
    fi
done