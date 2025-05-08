#!/bin/bash

set -e


YWDIR="/glade/derecho/scratch/jpan/b.e23.BMOM.f09_sx0.66av1.aqua.production.0530ywbranch"
CASENAME="b.e23.BMOM.f09_sx0.66av1.aqua.production.0621dlyocn"
CLONENAME="b.e23.BMOM.f09_sx0.66av1.aqua.production.0530ywbranch"
CASEDIR="/glade/u/home/jpan/work/MOM6_CASEDIRS/${CASENAME}"
CLONEDIR="/glade/u/home/jpan/work/MOM6_CASEDIRS/${CLONENAME}"
RUNDIR="/glade/derecho/scratch/jpan/${CASENAME}/run"
SRCROOT="/glade/u/home/jpan/work/CESM2.3_YOUWEI"
ARCHROOT="/glade/derecho/scratch/jpan/archive/"
RESTDATE="1014-10-01"

scp $SRCROOT/cime/config/cesm/machines/config_machines_try2b.xml $SRCROOT/cime/config/cesm/machines/config_machines.xml

$SRCROOT/cime/scripts/create_clone --case "${CASEDIR}" --clone "${CLONEDIR}"

cd "${CASEDIR}"

scp $CLONEDIR/user_nl_* ./

###No SourceMods

###./xmlchange NTASKS_ATM=512

###./xmlchange NTASKS_CPL=512

###./xmlchange NTASKS_LND=128

###./xmlchange NTASKS_ICE=384

###./xmlchange NTASKS_OCN=512

###./xmlchange ROOTPE_OCN=256

###./xmlchange ROOTPE_ICE=64

./xmlchange STOP_OPTION=nmonths

./xmlchange STOP_N=1

./xmlchange RESUBMIT=2

./xmlchange DEBUG=False


./xmlchange RUN_TYPE=branch

./xmlchange RUN_REFCASE="${CLONENAME}"

./xmlchange RUN_REFDATE="${RESTDATE}"

scp "${CLONEDIR}"/env_mach_pes.xml ./

./case.setup

scp "${ARCHROOT}"/"${CLONENAME}"/rest/"${RESTDATE}"-00000/* "${RUNDIR}"


###scp -r "${YWDIR}/run/INPUT/" "${RUNDIR}"
scp -r /glade/work/jpan/f09_sx0.66av1/ "${RUNDIR}"/INPUT

./xmlchange OCN_DOMAIN_MESH="${RUNDIR}"/INPUT/ESMF_mesh_sx0.66av1.nc

./xmlchange MASK_MESH="${RUNDIR}"/INPUT/ESMF_mesh_sx0.66av1.nc

./xmlchange CLM_FORCE_COLDSTART=on

./xmlchange ICE_DOMAIN_MESH="${RUNDIR}"/INPUT/ESMF_mesh_sx0.66av1.nc

./xmlchange LND_DOMAIN_FILE=domain.lnd.fv0.9x1.25_sx0.66av1.230330.nc

./xmlchange LND_DOMAIN_PATH=/glade/work/jpan/f09_sx0.66av1/

#qcmd -q main -A UPSU0063 -l walltime=00:20:00 ./case.build

./xmlchange --subgroup case.run PROJECT=UPSU0063

./xmlchange --subgroup case.st_archive PROJECT=UPSU0063

./xmlchange --subgroup case.run JOB_WALLCLOCK_TIME=2:00:00

./xmlchange --subgroup case.st_archive JOB_WALLCLOCK_TIME=1:00:00

#./case.submit
