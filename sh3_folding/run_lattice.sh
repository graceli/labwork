#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N foldit_max_q
#PBS -m e

# source gromacs just to be sure
. ~/.gmx

# stricter bash -- quits on error and unset variables
set -u
set -x

NPME=4
MDRUN=mdrun_d_mpi
GROMPP=grompp_d
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


function log {
	echo "$1: $2"
}

function run {
	MAXH=$1
	cd $base_dir/$PBS_ARRAYID
	echo "$PWD"
	cpt_file=""
	if [[ ! -e "prod.cpt" || ! -s "prod.cpt" ]]; then
		log "DEBUG" "checkpoint file does not exist or is empty .. starting a new run"
		# energy minimization	
		$GROMPP -f $PARAMS/em.mdp -c ${sysName} -p ${sysName} -o em.tpr
		$MDRUN -s em.tpr -deffnm em

		# equilibration
		$GROMPP -f $PARAMS/equil.mdp -c em.gro -p ${sysName}.top -o equil

		mpirun -np $NCORES $MDRUN -s equil.tpr -deffnm equil -npme $NPME
		
		# production runs
		$GROMPP -f $PARAMS/prod.mdp -c equil.gro -p ${sysName}.top -o prod.tpr
	else
		cpt_file="$base_dir/$PBS_ARRAYID/prod.cpt"
		log "INFO" "found checkpoint file $cpt_file"
	fi	

	log "DEBUG" "starting mdrun for $MAXH"
	mpirun -np $NCORES $MDRUN -s prod -deffnm prod -maxh $MAXH -cpt 720 -nosum -dlb auto -npme $NPME -cpo prod -cpi $cpt_file -noappend
}

# prepare environment
if [ "$PBS_ENVIRONMENT" != "PBS_INTERACTIVE" ]; then
  if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
  fi
fi

trap "exit 0" TERM KILL SIGINT SIGTERM EXIT

run $MAXH

