#!/bin/bash
#PBS -l nodes=1:compute-eth:ppn=8,walltime=48:00:00
#PBS -N foldit_mrsd

# stricter bash -- quits on error and unset variables
set -u
set -x
set -v
set -e

NPME=4
MDRUN=/project/p/pomes/cneale/GPC/exe/intel/gromacs-4.0.7/exec/bin/mdrun_openmpi-1.4.1
GROMPP=/project/p/pomes/cneale/GPC/exe/intel/gromacs-4.0.7/exec/bin/grompp
# NAME is passed in from the qsub 
sysName=${NAME}_${PBS_ARRAYID}
SHM=/dev/shm
base_dir=$PBS_O_WORKDIR
NRESUBMITS=0
#num=$NUM
DEBUG=0
MAXH=48
nodes=1
NCORES=16
PARAMS="../../params"

. ~/.gmx_407

function log {
	echo "INFO: $1"
}

function run {
	MAXH=$1
	cd $base_dir/$PBS_ARRAYID
	echo "$PWD"
	cpt_file=""
	if [[ ! -e "prod.cpt" || ! -s "prod.cpt" ]]; then
		echo "DEBUG: checkpoint file does not exist or is empty .. starting a new run"
		# energy minimization	
		$GROMPP -f $PARAMS/em.mdp -c ${sysName} -p ${sysName} -o em.tpr
		$MDRUN -s em.tpr -deffnm em

		# equilibration
		log "starting equilibration stage"
		$GROMPP -f $PARAMS/equil_npt.mdp -c em.gro -p ${sysName}.top -o equil_npt

		mpirun -np $NCORES $MDRUN -s equil_npt.tpr -deffnm equil_npt -nosum -dlb auto -npme $NPME
		#$MDRUN -s equil.tpr -deffnm equil

	
		log "making tpr for production run"
		# production runs
		$GROMPP -f $PARAMS/prod_npt.mdp -c equil_npt.gro -p ${sysName}.top -o prod_npt
	else
		cpt_file="$base_dir/$PBS_ARRAYID/prod.cpt"
		echo "INFO: found checkpoint file $cpt_file"
	fi	

	echo "DEBUG: starting mdrun for $MAXH"
	mpirun -np $NCORES $MDRUN -s prod_npt -deffnm prod_npt -maxh $MAXH -cpt 720 -nosum -dlb auto -npme $NPME -cpo prod -cpi $cpt_file -noappend
}

# prepare environment
if [ "$PBS_ENVIRONMENT" != "PBS_INTERACTIVE" ]; then
  if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
  fi
fi

#trap "exit 0" TERM KILL SIGINT SIGTERM EXIT

run $MAXH

# if [ "$DEBUG" -eq "1" ]; then
#    DEBUG_FLAGS="-l nodes=1:ppn=8,walltime=02:00:00 -q debug"
# else
#    DEBUG_FLAGS=""
# fi
