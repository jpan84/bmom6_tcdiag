#!/bin/bash -l
#PBS -N ncremap_hy2pres
#PBS -A UCIS0005
#PBS -l walltime=00:55:00
#PBS -l select=2:ncpus=8:mem=180GB
#PBS -q casper

set -e
module load nco

###EXP_DIR = "$1"
cd "${EXP_DIR}"
pwd

export diri="./hist/"
export diro="./hist_onpres_gnupar/"
tapes=("h0a")
export vrt_in="/glade/u/home/jpan/grids/L26_hyb.nc"
export vrt_out="/glade/u/home/jpan/grids/L26_plevs.nc"

if [ ! -d "$diro" ]; then
  mkdir -p "$diro"
fi

for tp in "${tapes[@]}"; do
  # Find files and pipe them to parallel
  find "${diri}" -name "*.${tp}.*.nc" | parallel -j 16 \
    "ncremap --no_stdin --thr_nbr=1 --vrb=2 --vrt_in=${vrt_in} --vrt_out=${vrt_out} --ps_nm={}/PS {} ${diro}/{/.}.onpres.nc"

done
