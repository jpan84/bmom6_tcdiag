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
  files=("${dir}"/*.nc)
  if [[ ! -e "${files[0]}" ]]; then
    continue
  fi
  for nc in "${files[@]}"; do
    #xz "${nc}"
    ###xz too expensive to run on interactive node
    qcmd -q main -l walltime=00:02:00 -A UPSU0063 xz -v "${nc}"
  done
done
