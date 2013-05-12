#!/bin/sh

#PBS -A uix-840-af
#PBS -N extend
#PBS -l walltime=48:00:00
#PBS -l nodes=7:ppn=8
#PBS -l flags=restartable
#PBS -l gattr=ckpt
#PBS -o MyExtend-%I.out
#PBS -e MyExtend-%I.err

module load compilers/intel/11.1.059
module load mpi/openmpi/1.3.4_intel
echo "Got $NSLOTS processors."
. /home/grace/.gmx

#set -u
set -e
set -x

NPME=-1
sysName=sys${MOAB_JOBARRAYINDEX}
NRESUBMITS=2
num=$NUM
base_dir=$PWD

function run {
	MAXH=$1
	cpt_file=sys${MOAB_JOBARRAYINDEX}_prod.cpt
	cd $base_dir/${MOAB_JOBARRAYINDEX}
	echo 'DEBUG: starting mdrun for $MAXH'

	mpirun mdrun -s ${sysName}_prod -deffnm ${sysName}_prod -maxh $MAXH -cpt 720 -nosum -dlb auto -npme $NPME -cpo ${sysName}_prod -cpi $cpt_file
}

run 46

if [ "$num" -lt "$NRESUBMITS" ]; then
	num=$((NUM+1))
	echo "resubmitting - sequence $num for replica $MOAB_JOBARRAYINDEX"
	qsub -v NUM=$num,MOAB_JOBARRAYINDEX=${MOAB_JOBARRAYINDEX} -N ${JOB_NAME}_${MOAB_JOBARRAYINDEX}_${num} /home/grace/labwork/abeta_analysis/scripts/rerun_colosse.sh
fi
