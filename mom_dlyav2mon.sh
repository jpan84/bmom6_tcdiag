#!/bin/bash

set -e

module load nco

CASENAME="b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.250121_unseed_all"
tapes=("hm" "h")

for tape in "${tapes[@]}"; do
    echo $tape
    MOSTRECENT=`ls --color=never ./$CASENAME.mom6.$tape.*.nc -r | head -n 1`
    ###echo $MOSTRECENT
    LASTDATE=$(echo "$MOSTRECENT" | awk -F. '{print $(NF-1)}')
    ###echo $LASTDATE
    IFS='-' read -r LASTYYYY LASTMM LASTDD <<< "$LASTDATE"
    echo $LASTDD

    ###noleap
    DDINMO=$(cal $LASTMM 0001 | awk 'NF {DAYS = $NF}; END {print DAYS}')
    ###echo $DDINMO

    if [[ "$LASTDD" == "$DDINMO" ]]; then
        echo "Month ${LASTYYYY}-${LASTMM} is complete. Averaging and archiving tape ${tape}." 
        MOFILES="$CASENAME.mom6.$tape.$LASTYYYY-$LASTMM-*.nc"
        ncra ${MOFILES} "$CASENAME.mom6.${tape}cust_avg.$LASTYYYY-$LASTMM.nc"
        #mv -v "*cust_avg.nc" "${ARCHROOT}/${CASENAME}/ocn/hist"
        rm -v $MOFILES
    else
        echo "Month ${LASTYYYY}-${LASTMM} is incomplete. Not averaging."
    fi
done

