###process high-freq 2D atm using NodeFileFilter masks and yhourmean output
#!/bin/bash -l

set -e
module load nco
module load cdo

ARCHDIR="/glade/campaign/univ/upsu0032/jpan_tcfields"
CASES=("b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/" "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/" "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/")

DATDIR="hist_0010_h1i/"
DATFIL="*.h1i.*.nc"
MSKDIR=$DATDIR
MSKEXT=".nff_2mps"
MSKNAM="TC_R2"

NORMF="hist_norms/yhourmean_h1i.nc"
DELVARS="ch4vmr,co2vmr,date,datesec,f11vmr,f12vmr,n2ovmr,ndcur,nscur,nsteph,sol_tsi"

cd $ARCHDIR

if [ "$(pwd)" != "$ARCHDIR" ]; then
  echo "Not running in the correct directory. Exiting to prevent destroying things."
  exit 1
fi



for case in "${CASES[@]}"; do
  cd $ARCHDIR/$case
  mkdir -p ${MSKNAM}_masked
  mkdir -p yhoureddy
  mkdir -p yhoureddy_${MSKNAM}_masked

  DATS=$(find "${DATDIR}" -maxdepth 1 -name "${DATFIL}" | sort)
  
  for datf in $DATS; do
    echo "Working on ${datf}"
    cdo mul -delname,${DELVARS} ${datf} -selname,${MSKNAM} ${MSKDIR}/$(basename ${datf})${MSKEXT} "${MSKNAM}_masked/$(basename ${datf})"
    ncks -O -4 -L 1 "${MSKNAM}_masked/$(basename ${datf})" "${MSKNAM}_masked/$(basename ${datf})"
  done

done

exit
##################################################################

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
