#!/bin/bash -l

set -e
module load nco

cd /glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/atm
pwd

diri="./hist_regrid_0.25x0.25/"
diro="./hist_regrid_0.25x0.25_onpres/"
tapes=("h0a")
vrt_in="/glade/u/home/jpan/grids/L26_hyb.nc"
vrt_out="/glade/u/home/jpan/grids/L26_plevs.nc"

if [ ! -d "$diro" ]; then
  mkdir -p "$diro"
fi

for tp in "${tapes[@]}"; do
  for nc in "${diri}"/*."${tp}".*.nc; do
    echo $nc
    bn=$(basename "${nc}")
    ncremap --no_stdin --thr_nbr=1 --vrb=2 --vrt_in="${vrt_in}" --vrt_out="${vrt_out}" --ps_nm="${nc}"/PS "${nc}" "${diro}"/"${bn}"
  done
done
