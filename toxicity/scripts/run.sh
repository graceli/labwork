#!/bin/bash
#PBS -q qwork@mp2 
#PBS -l walltime=120:00:00 -l nodes=2:ppn=1
#PBS -N tox_ino

# stricter bash -- quits on error and unset variables
set -u
set -e

NPME=-1
#EXE=/home/nealechr/exe/gromacs-4.5.4/exec/bin/mdrun_mpi
EXE=mdrun_mpi
sysName=system_${PBS_ARRAYID}_final
SHM=$SCRATCH
base_dir=$PBS_O_WORKDIR
NRESUBMITS=0
#num=$NUM
DEBUG=0
MAXH=100
nodes=2
NCORES=$(echo "$nodes*24" | bc)
PARAMS="../params"


function log {
	echo "INFO: $1"
}

function clean_exit {
	echo "DEBUG: starting clean_exit ... on ARRAYID $PBS_ARRAYID"
	cd $SHM
	cp -p *.log *.tpr *.xtc *.edr *.gro *.cpt $base_dir/$PBS_ARRAYID
	cd $base_dir
	rm -rf $SHM
	echo "DEBUG: cleaned up $SHM"
}

function run {
	MAXH=$1
	cd $base_dir/$PBS_ARRAYID 
	echo "$PWD"
	cpt_file=""
	if [[ ! -e "prod.cpt" || ! -s "prod.cpt" ]]; then
		echo "DEBUG: checkpoint file does not exist or is empty .. starting a new run"
		#tmpname=`mktemp -d backup.XXX`
		#mv *.cpt $tmpname
		
		# workflow	
		# energy minimization	
		echo "grompp -f $PARAMS/em.mdp -c ${sysName} -p ${sysName} -o em.tpr"	
		grompp_mpi -f $PARAMS/em.mdp -c ${sysName} -p ${sysName} -o em.tpr
		echo "after grompp"
		echo "mdrun_mpi -s em.tpr -deffnm em"
		mdrun_mpi -s em.tpr -deffnm em
		echo "after em"
		# equilibration
		grompp_mpi -f $PARAMS/equil.mdp -c em.gro -p ${sysName}.top -o equil
		echo "mpirun -np $NCORES $EXE -s equil.tpr -deffnm equil -npme -1"
		mpirun -np $NCORES $EXE -s equil.tpr -deffnm equil -npme -1 
		
		# production runs
		grompp_mpi -f $PARAMS/prod.mdp -c equil.gro -p ${sysName}.top -o prod.tpr
	else
		cpt_file="$base_dir/$PBS_ARRAYID/prod.cpt"
		echo "INFO: found checkpoint file $cpt_file"
	fi	

	echo "DEBUG: starting mdrun for $MAXH"
	mpirun -np $NCORES $EXE -s prod -deffnm prod -maxh $MAXH -cpt 720 -nosum -dlb auto -npme $NPME -cpo prod -cpi $cpt_file -noappend
}

#if [ "$RET" -gt "0" ]; then
#	`date` >>ERROR
#	echo "there was a problem" >> ERROR
#	exit
#fi

# prepare environment
if [ "$PBS_ENVIRONMENT" != "PBS_INTERACTIVE" ]; then
  if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
  fi
fi

trap "clean_exit; exit 0" TERM KILL SIGINT SIGTERM EXIT

run $MAXH

if [ "$DEBUG" -eq "1" ]; then
    DEBUG_FLAGS="-l nodes=1:ppn=8,walltime=02:00:00 -q debug"
else
    DEBUG_FLAGS=""
fi

#clean_exit
