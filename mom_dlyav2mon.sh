#!/bin/bash -l

set -e

module load nco

ARCHROOT="/glade/derecho/scratch/jpan/archive/"
CASENAME="REPLCASE"
tapes=("hm" "h")

cd "$ARCHROOT/$CASENAME/ocn/hist"

for tape in "${tapes[@]}"; do
    echo $tape
    MOSTRECENT=`ls --color=never ./$CASENAME.mom6.$tape.*.nc -r | head -n 1`
    ###echo $MOSTRECENT
    LASTDATE=$(echo "$MOSTRECENT" | awk -F. '{print $(NF-1)}')
    ###echo $LASTDATE
    IFS='-' read -r LASTYYYY LASTMM LASTDD <<< "$LASTDATE"
    ###echo $LASTDD

    while true; do
        ###if no diag files exist for the month, then break
        MOFILES="$CASENAME.mom6.$tape.$LASTYYYY-$LASTMM-*.nc"
        ###ls -1 ${MOFILES}
        if [[ -z "$(ls -1 ${MOFILES} 2>/dev/null)" ]]; then
            break
        fi

        DDCOUNT=$(ls -1 ${MOFILES} 2>/dev/null | wc -l)

        ###noleap
        DDINMO=$(cal $LASTMM 0001 | awk 'NF {DAYS = $NF}; END {print DAYS}')
        ###echo $DDINMO
    
        if [[ "$DDCOUNT" == "$DDINMO" ]]; then
            echo "Month ${LASTYYYY}-${LASTMM} is complete. Averaging and archiving tape ${tape}." 
            ncra ${MOFILES} "$CASENAME.mom6.${tape}cust_avg.$LASTYYYY-$LASTMM.nc"
            #mv -v "*cust_avg.nc" "${ARCHROOT}/${CASENAME}/ocn/hist"
            rm -v $MOFILES
        else
            echo "Month ${LASTYYYY}-${LASTMM} is incomplete. Not averaging."
        fi

        ###convert back to ints
        LASTYYYY=$((10#$LASTYYYY))
        LASTMM=$((10#$LASTMM))
        ###echo $LASTMM

        ###decrement month and retain leading zeroes
        if [[ $LASTMM -eq 1 ]]; then
            ((LASTYYYY--))
            LASTMM=12 
        else
            ((LASTMM--))            
        fi
        LASTMM=$(printf "%02d" "$LASTMM")
        LASTYYYY=$(printf "%04d" "$LASTYYYY")

    done
done

