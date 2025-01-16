#!/bin/bash

set -e

CASENAME="b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.250116_unseed_all"
CLONENAME="b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.1121_ctrl"
CASEDIR="/glade/u/home/jpan/work/MOM6_CASEDIRS/${CASENAME}"
CLONEDIR="/glade/u/home/jpan/work/MOM6_CASEDIRS/${CLONENAME}"
RUNDIR="/glade/derecho/scratch/jpan/${CASENAME}/run"
SRCROOT="/glade/u/home/jpan/work/CESM23a17d_MOM/"
ARCHROOT="/glade/derecho/scratch/jpan/archive/"
RESTDATE="0001-02-01"

SEEDDRVDIR="/glade/u/home/jpan/work/betacast-pivot/py_atm_to_cam/tcseed/"

$SRCROOT/cime/scripts/create_clone --case "${CASEDIR}" --clone "${CLONEDIR}"

cd "${CASEDIR}"

###Disable restart checksum
scp /glade/u/home/jpan/work/MOM6_CASEDIRS/b.e23.BMOM.ne120pg3_sx0.66av1.aqua.zmtau.0826/user_nl_mom ./

scp /glade/u/home/jpan/aquaptc/mom6diag/1201_diag_table_alldaily ./SourceMods/src.mom/diag_table

sed -i "s/CASENAME/$CASENAME/g" ./SourceMods/src.mom/diag_table

./xmlchange STOP_OPTION=ndays

./xmlchange STOP_N=1

./xmlchange RESUBMIT=0

./xmlchange DEBUG=False

./xmlchange RUN_TYPE=branch

./xmlchange RUN_REFCASE="${CLONENAME}"

./xmlchange RUN_REFDATE="${RESTDATE}"

mkdir -p $RUNDIR

scp "${ARCHROOT}"/"${CLONENAME}"/rest/"${RESTDATE}"-00000/* "${RUNDIR}"

ocnrest="/glade/derecho/scratch/jpan/archive/b.e23.BMOM.f09_sx0.66av1.aqua.production.0930_a17d_branch/rest/0005-02-01-00000/b.e23.BMOM.f09_sx0.66av1.aqua.production.0930_a17d_branch.mom6.r.0005-02-01-00000.nc"

scp ${ocnrest} ${RUNDIR}/${CLONENAME}.mom6.r.${RESTDATE}-00000.nc

scp -r /glade/work/zarzycki/CESM_files/mom/INPUT/ "${RUNDIR}"/INPUT

cat > user_nl_cam <<EOF
nhtfrq=0,-6,-1
mfilt=0,1,24
fincl2='UBOT:I','VBOT:I','Z300:I','SST:I','Q850:I','PS:I','PSDRY:I','PRECT:I','U850:I','V850:I','PSL:I','TBOT:I','T850:I','UBOT:I','VBOT:I','TMQ:I','Z500:I','U200:I','V200:I','OMEGA500:I','OMEGA850:I','T200:    I','T500:I','U500:I','V500:I', 'TAUX:I', 'TAUY:I'
fincl3='PRECT:A','PRECC:A'
se_write_restart_unstruct              = .true.
ncdata='/glade/campaign/cesm/cesmdata/inputdata/atm/cam/inic/se/ape_cam4_ne30np4_L26_c170417.nc'
prescribed_strataero_feedback             = .false.
ch4vmr            = 1.650e-6
n2ovmr    = 0.306e-6
co2vmr    = 348.0e-6
solar_const = 1365.0
omega = 7.292115e-5
prescribed_ozone_file     = 'apeozone_cam3_5_54.nc'
prescribed_ozone_datapath = '/glade/campaign/cesm/cesmdata/inputdata/atm/cam/ozone'
prescribed_ozone_name     = 'OZONE'
prescribed_ozone_type     = 'CYCLICAL'
prescribed_ozone_cycle_yr = 1990
prescribed_aero_file = ''
sday   = 86164.10063718943
rearth = 6.37100e6
gravit = 9.79764
mwdry  = 28.96623324623746
mwh2o  = 18.01618112892741
cpwv   = 1.846e3
rad_climate = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2', 'N:ozone:O3', 'N:N2O:N2O', 'N:CH4:CH4', 'N:CFC11:CFC11','N:CFC12:CFC12'
use_topo_file             = .false.
aerodep_flx_file               = ''
scale_dry_air_mass             =  101080.0D0
atm_dep_flux           = .false.
EOF

cat > user_nl_clm <<EOF
fsurdat = '/glade/work/zarzycki/CESM_files/fsurdat/surfdata_BAQUA_CLM5_ne30np4.c241120.nc'
urban_traffic = .false.
EOF

cat <<EOF > config.unseed
DRIVER_DIRECTORY = $SEEDDRVDIR
path_to_case = $CASEDIR
path_to_rundir = $RUNDIR
MAXTIME=00310201
TEMPESTEXTREMESDIR = /glade/work/zarzycki/derecho/tempestextremes_noMPI
CONNECTDAT = "/glade/u/home/zarzycki/tempest-scripts/hyperion/ne30.connect_v2.dat"
TEMPESTFILE = $RUNDIR/cyclones_tempest
TEMPESTTMP = $SEEDDRVDIR/cyclones_tempest.tmp
HTRACKSTR = "h1i"
vortex_namelist = /glade/u/home/jpan/work/betacast-pivot/py_atm_to_cam/tcseed/unseed.ape.nl.100.bmom
do_seed = false
st_archive_ontheway = true
KILLSTR = "UNSEED_DONE"
EOF

scp config.unseed config.unseed.backup

./xmlchange JOB_PRIORITY=economy

./xmlchange --subgroup case.run PROJECT=UPSU0032

./xmlchange --subgroup case.st_archive PROJECT=UPSU0032

./xmlchange --subgroup case.run JOB_WALLCLOCK_TIME=00:15:00

./xmlchange --subgroup case.st_archive JOB_WALLCLOCK_TIME=00:15:00

./case.build --clean

./case.setup

qcmd -q main -A UPSU0032 -l walltime=00:30:00 ./case.build

cd $SEEDDRVDIR
###scp tc-seed-driver.sh.ne30.1121 tc-seed-driver.sh

###delete daily restarts from archive
scp /glade/u/home/jpan/aquaptc/mom6diag/del_st_rest.sh $RUNDIR
sed -i "s/CASENAME/$CASENAME/g" $RUNDIR/del_st_rest.sh
(crontab -l 2>/dev/null; echo "0 0,12 * * * ${RUNDIR}/del_st_rest.sh") | crontab -

###./case.submit
./tc-seed-driver.sh.ne30.1121 $CASEDIR/config.unseed
