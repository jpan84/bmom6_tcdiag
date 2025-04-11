###compress archived restarts
###Run on derecho
#!/bin/bash -l

set -e

###no trailing slash
ARCHREST="/glade/derecho/scratch/jpan/archive/REPLCASE/rest"

cd $ARCHREST

if [ "$(pwd)" != "$ARCHREST" ]; then
  echo "Not running in the correct directory. Exiting to prevent destroying things."
  exit 1
fi

for dir in "${ARCHREST}"/*/; do
  ###rm "${dir}"/*.cam.h*.nc || true
  ###rm "${dir}"/*.clm2.h0*.nc || true
  for nc in "${dir}"/*.nc; do
    #xz "${nc}"
    ###xz too expensive to run on interactive node
    qcmd -q main -l walltime=00:01:00 -A UPSU0063 xz -v "${nc}"
    ###if [[ $? -eq 0 ]]; then # check if xz was successful
    ###  rm "${nc}"
    ###fi
  done
done
