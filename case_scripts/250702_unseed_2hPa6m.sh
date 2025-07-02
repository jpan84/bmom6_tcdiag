###Production run of aggressive unseeding in ne120np4
#!/bin/bash

set -e

CASENAME="b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250702_unseed2hPa6m"
CLONENAME="b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250206_1degeqm885"
CASEDIR="/glade/u/home/jpan/work/MOM6_CASEDIRS/${CASENAME}"
CLONEDIR="/glade/u/home/jpan/work/MOM6_CASEDIRS/${CLONENAME}"
RUNDIR="/glade/derecho/scratch/jpan/${CASENAME}/run"
SRCROOT="/glade/u/home/jpan/work/CESM23a17d_MOM/"
ARCHROOT="/glade/derecho/scratch/jpan/archive/"
DIAGSRC="/glade/u/home/jpan/aquaptc/bmom6_tcdiag/"
SEEDDRVDIR="/glade/u/home/jpan/work/betacast-pivot/py_atm_to_cam/tcseed/"
STARTMO="02"
RESTDATE="0005-${STARTMO}-01"
USEPROJ="UPSU0063"

$SRCROOT/cime/scripts/create_clone --case "${CASEDIR}" --clone "${CLONEDIR}"

cd "${CASEDIR}"

scp ${DIAGSRC}/250415_diag_table_dlyrest ./SourceMods/src.mom/diag_table

sed -i "s/CASENAME/$CASENAME/g" ./SourceMods/src.mom/diag_table

#add Z850
scp ${DIAGSRC}/cam_diagnostics.F90 ./SourceMods/src.cam/

#Disable ocean state checksum and include rho2 diagnostics
scp ${DIAGSRC}/user_nl_mom ./

#Make sure mom6 hm, rho2, sfc files are moved to archive
scp ${DIAGSRC}/env_archive.xml ./

scp ${DIAGSRC}/pkg_cldoptics.F90.bacmeister ./SourceMods/src.cam/pkg_cldoptics.F90

./xmlchange STOP_OPTION=ndays

./xmlchange STOP_N=1

./xmlchange RESUBMIT=0

./xmlchange DEBUG=False

./xmlchange RUN_TYPE=branch

./xmlchange RUN_REFCASE="${CLONENAME}"

./xmlchange RUN_REFDATE="${RESTDATE}"

mkdir -p $RUNDIR

scp "${ARCHROOT}"/"${CLONENAME}"/rest/"${RESTDATE}"-00000/* "${RUNDIR}"

###ocnrest="/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.250115_h79l90/rest/0061-02-01-00000/b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.250115_h79l90.mom6.r.0061-02-01-00000.nc"

###scp ${ocnrest} ${RUNDIR}/${CLONENAME}.mom6.r.${RESTDATE}-00000.nc

scp -r /glade/work/zarzycki/CESM_files/mom/INPUT/ "${RUNDIR}"/INPUT
scp /glade/u/home/jpan/ocean_rho2_250415.nc "${RUNDIR}"/INPUT

cat > user_nl_cam <<EOF
nhtfrq=0,-6,-1
mfilt=0,1,24
fincl1='CAPE:A','VZ:A','UWzm:A'
fincl2='UBOT:I','VBOT:I','Z300:I','SST:I','Q850:I','Q500:I','Q200:I','PS:I','PSDRY:I','PRECT:I','U850:I','V850:I','PSL:I','TBOT:I','T850:I','UBOT:I','VBOT:I','TMQ:I','Z500:I','U200:I','V200:I','OMEGA500:I','OMEGA850:I','T200:I','T500:I','U500:I','V500:I', 'TAUX:I', 'TAUY:I', 'FLUT:I', 'Z850:I', 'CAPE:I'
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
./xmlchange NTASKS_ATM=2816
./xmlchange NTHRDS_ATM=1
./xmlchange ROOTPE_ATM=0

./xmlchange NTASKS_CPL=2816
./xmlchange NTHRDS_CPL=1
./xmlchange ROOTPE_CPL=0

./xmlchange NTASKS_GLC=128
./xmlchange NTHRDS_GLC=1
./xmlchange ROOTPE_GLC=2816

./xmlchange NTASKS_ROF=1920
./xmlchange NTHRDS_ROF=1
./xmlchange ROOTPE_ROF=128

./xmlchange NTASKS_ICE=2688
./xmlchange NTHRDS_ICE=1
./xmlchange ROOTPE_ICE=128

./xmlchange NTASKS_LND=128
./xmlchange NTHRDS_LND=1
./xmlchange ROOTPE_LND=0

./xmlchange NTASKS_OCN=128
./xmlchange NTHRDS_OCN=1
./xmlchange ROOTPE_OCN=2816


cat > config.unseed << EOF
DRIVER_DIRECTORY = $SEEDDRVDIR
path_to_case = $CASEDIR
path_to_rundir = $RUNDIR
MAXTIME=00060201
TEMPESTEXTREMESDIR = /glade/work/zarzycki/derecho/tempestextremes_noMPI
CONNECTDAT = "/glade/u/home/jpan/ne120np4_connect_v2.dat"
TEMPESTFILE = $RUNDIR/cyclones_tempest
TEMPESTTMP = $SEEDDRVDIR/cyclones_tempest.tmp
HTRACKSTR = "h1i"
vortex_namelist = /glade/u/home/jpan/work/betacast-pivot/py_atm_to_cam/tcseed/unseed.ape.nl.025.1
do_seed = false
st_archive_ontheway = true
KILLSTR = "UNSEED_DONE"
EOF

./case.setup --reset
./case.setup
./case.build --clean-all

./xmlchange JOB_PRIORITY=economy

./xmlchange --subgroup case.run PROJECT=$USEPROJ

./xmlchange --subgroup case.st_archive PROJECT=$USEPROJ

./xmlchange --subgroup case.run JOB_WALLCLOCK_TIME=0:07:00

./xmlchange --subgroup case.st_archive JOB_WALLCLOCK_TIME=0:05:00

qcmd -q main -A $USEPROJ -l walltime=00:30:00 ./case.build

###create calendar month averages from daily averages and delete daily averages
###make sure output tapes in the script match diag_table
scp ${DIAGSRC}/mom_dlyav2mon.sh $RUNDIR
sed -i "s/REPLCASE/$CASENAME/g" $RUNDIR/mom_dlyav2mon.sh
ssh -n cron "(crontab -l 2>/dev/null; echo \"0 */3 * * * ssh -n casper '/bin/bash ${RUNDIR}/mom_dlyav2mon.sh &> ${RUNDIR}/mom_dlyav2mon_cron.out'\") | crontab -"

###apply ncks lossless compression to history archive
scp ${DIAGSRC}/ncks_hist.sh $RUNDIR
sed -i "s/REPLCASE/$CASENAME/g" $RUNDIR/ncks_hist.sh
ssh -n cron "(crontab -l 2>/dev/null; echo \"45 */3 * * * ssh -n casper '/bin/bash ${RUNDIR}/ncks_hist.sh &> ${RUNDIR}/ncks_hist_cron.out'\") | crontab -"

######compress short-term restarts with xz
###scp ${DIAGSRC}/xz_st_rest.sh $RUNDIR
###sed -i "s/REPLCASE/$CASENAME/g" $RUNDIR/xz_st_rest.sh
###ssh -n cron "(crontab -l 2>/dev/null; echo \"0 */6 * * * ssh -n derecho '/bin/bash ${RUNDIR}/xz_st_rest.sh &> ${RUNDIR}/xz_st_rest_cron.out'\") | crontab -"

scp ${DIAGSRC}/del_st_rest.sh $RUNDIR
sed -i "s/REPLCASE/$CASENAME/g" $RUNDIR/del_st_rest.sh
sed -i "s/REPLMO/$STARTMO/g" $RUNDIR/del_st_rest.sh
ssh -n cron "(crontab -l 2>/dev/null; echo \"0 */3 * * * ssh -n derecho '/bin/bash ${RUNDIR}/del_st_rest.sh &> ${RUNDIR}/del_st_rest_cron.out'\") | crontab -"

###scp ${DIAGSRC}/del_st_rest.sh $RUNDIR
###sed -i "s/REPLCASE/$CASENAME/g" $RUNDIR/del_st_rest.sh
###sed -i "s/REPLMO/$STARTMO/g" $RUNDIR/del_st_rest.sh
###(crontab -l 2>/dev/null; echo "0 0,12 * * * ${RUNDIR}/del_st_rest.sh") | crontab -

echo "Cron jobs created to compute MOM monthly averages, compress history, and delete st restarts."
echo "Run crontab -e to view or remove these jobs."

cd $SEEDDRVDIR
###./case.submit
./tc-seed-driver.sh.250702_unseed_sznl $CASEDIR/config.unseed
