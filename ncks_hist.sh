###compress (lossless) hist by ~1/4 to 1/3
###by rechunking with ncks
#!/bin/bash

set -e
module load nco

#Omit trailing slash to allow pwd check to work
ARCHDIR="/glade/derecho/scratch/jpan/archive/REPLCASE"
comps=("atm/" "lnd/" "ocn/" "ice/")

cd $ARCHDIR

if [ "$(pwd)" != "$ARCHDIR" ]; then
  echo "Not running in the correct directory. Exiting to prevent destroying things."
  exit 1
fi

for cmp in "${comps[@]}"; do
  ###echo "${ARCHDIR}/${cmp}/hist/*.nc"
  for nc in "${ARCHDIR}"/"${cmp}"/hist/*.nc; do
    if [[ "${nc}" == *".mom6.h."* || "${nc}" == *".mom6.hm."* || "${nc}" == *".mom6.rho2."* ]]; then
      continue
    fi
    if ncdump -h ${nc} | grep ":history" | grep -q "ncks -O -4 -L 1"; then
      ###echo "Skipping ${nc}"
      continue
    else
      ###echo "Rechunking ${nc}"
      ncks -O -4 -L 1 "${nc}" "${nc}"
    fi
  done
done
