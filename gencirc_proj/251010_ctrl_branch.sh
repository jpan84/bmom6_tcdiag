#!/bin/bash

set -e

CASENAME="b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr"
CLONENAME="b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl"
CASEDIR="/glade/u/home/jpan/work/MOM6_CASEDIRS/${CASENAME}"
CLONEDIR="/glade/u/home/jpan/work/MOM6_CASEDIRS/${CLONENAME}"
RUNDIR="/glade/derecho/scratch/jpan/${CASENAME}/run"
SRCROOT="/glade/u/home/jpan/work/CESM23a17d_MOM/"
ARCHROOT="/glade/campaign/univ/upsu0032/jpan_aquaptc/"
DIAGSRC="/glade/u/home/jpan/aquaptc/bmom6_tcdiag/"
STARTMO="02"
RESTDATE="0012-${STARTMO}-01"
USEPROJ="UPSU0032"

$SRCROOT/cime/scripts/create_clone --case "${CASEDIR}" --clone "${CLONEDIR}"

cd "${CASEDIR}"

scp ${DIAGSRC}/250415_diag_table_monrest ./SourceMods/src.mom/diag_table

sed -i "s/CASENAME/$CASENAME/g" ./SourceMods/src.mom/diag_table

#add Z850
scp ${DIAGSRC}/cam_diagnostics.F90 ./SourceMods/src.cam/

#Disable ocean state checksum and include rho2 diagnostics
scp ${DIAGSRC}/user_nl_mom ./

#Make sure mom6 hm, rho2, sfc files are moved to archive
scp ${DIAGSRC}/env_archive.xml ./

scp ${DIAGSRC}/pkg_cldoptics.F90.bacmeister ./SourceMods/src.cam/pkg_cldoptics.F90

./xmlchange STOP_OPTION=nmonths

./xmlchange STOP_N=1

./xmlchange RESUBMIT=23

./xmlchange DEBUG=False

./xmlchange RUN_TYPE=branch

./xmlchange RUN_REFCASE="${CLONENAME}"

./xmlchange RUN_REFDATE="${RESTDATE}"

mkdir -p $RUNDIR

scp "${ARCHROOT}"/"${CLONENAME}"/rest/"${RESTDATE}"-00000/* "${RUNDIR}"

###ocnrest="/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.250130_h80l89/rest/0051-02-01-00000/b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.250130_h80l89.mom6.r.0051-02-01-00000.nc"

###scp ${ocnrest} ${RUNDIR}/${CLONENAME}.mom6.r.${RESTDATE}-00000.nc

scp -r /glade/work/zarzycki/CESM_files/mom/INPUT/ "${RUNDIR}"/INPUT
scp /glade/u/home/jpan/ocean_rho2_250415.nc "${RUNDIR}"/INPUT

cat > user_nl_cam <<EOF
nhtfrq=0,-6,-1
mfilt=0,1,24
fincl1='CAPE:A'
fincl2='Z200:I','Z300:I','SST:I','Q850:I','Q500:I','PS:I','PRECT:I','U850:I','V850:I','T850:I','UBOT:I','VBOT:I','TMQ:I','Z500:I','U200:I','V200:I','OMEGA200:I','OMEGA500:I','OMEGA850:I','T200:I','T500:I','U500:I','V500:I', 'TAUX:I', 'TAUY:I', 'FLUT:I', 'Z850:I', 'CAPE:I','IVT:I','VU:I','VV:I','VT:I','VQ:I','VZ:I','UWzm:I','OMEGAT:I','U:I','V:I','OMEGA:I','T:I','Q:I','Z3:I'
fincl3='PRECT:A','PRECC:A'
se_write_restart_unstruct              = .true.
ncdata='/glade/work/zarzycki/CESM_files/ncdata/ncdata_ne120np4_L26_aqua.nc'
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
zmconv_tau = 1800.0
cldfrc_rhminh          =  0.800D0
cldfrc_rhminl          =  0.885D0
EOF

#change processor layout
./xmlchange NTASKS_ATM=2944
./xmlchange NTHRDS_ATM=1
./xmlchange ROOTPE_ATM=0

./xmlchange NTASKS_CPL=768
./xmlchange NTHRDS_CPL=1
./xmlchange ROOTPE_CPL=128

./xmlchange NTASKS_ICE=2048
./xmlchange NTHRDS_ICE=1
./xmlchange ROOTPE_ICE=896

./xmlchange NTASKS_LND=128
./xmlchange NTHRDS_LND=1
./xmlchange ROOTPE_LND=0

./xmlchange NTASKS_OCN=128
./xmlchange NTHRDS_OCN=1
./xmlchange ROOTPE_OCN=2944

./case.setup --reset
./case.setup
./case.build --clean-all

./xmlchange JOB_PRIORITY=economy

./xmlchange --subgroup case.run PROJECT=$USEPROJ

./xmlchange --subgroup case.st_archive PROJECT=$USEPROJ

./xmlchange --subgroup case.run JOB_WALLCLOCK_TIME=0:40:00

./xmlchange --subgroup case.st_archive JOB_WALLCLOCK_TIME=0:10:00

qcmd -q main -A $USEPROJ -l walltime=00:30:00 ./case.build

./case.submit

###apply ncks lossless compression to history archive
scp ${DIAGSRC}/ncks_hist.sh $RUNDIR
sed -i "s/REPLCASE/$CASENAME/g" $RUNDIR/ncks_hist.sh
ssh -n cron "(crontab -l 2>/dev/null; echo \"45 */3 * * * ssh -n casper '/bin/bash ${RUNDIR}/ncks_hist.sh &> ${RUNDIR}/ncks_hist_cron.out'\") | crontab -"

scp ${DIAGSRC}/del_st_rest.sh $RUNDIR
sed -i "s/REPLCASE/$CASENAME/g" $RUNDIR/del_st_rest.sh
sed -i "s/REPLMO/$STARTMO/g" $RUNDIR/del_st_rest.sh
ssh -n cron "(crontab -l 2>/dev/null; echo \"0 */3 * * * ssh -n casper '/bin/bash ${RUNDIR}/del_st_rest.sh &> ${RUNDIR}/del_st_rest_cron.out'\") | crontab -"
echo "Cron job created to compress history and delete st restarts."
echo "Run crontab -e on cron server to view or remove these jobs."
