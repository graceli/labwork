#!/bin/sh
#$ -N test
#$ -P uix-840-ac
#$ -A uix-840-ac
##$ -l h_rt=48:00:00
##$ -pe default 48 
##$ -q med
#$ -S /bin/bash
#$ -cwd
#$ -notify

module load compilers/intel/11.1.059
module load mpi/openmpi/1.3.4_intel
export OMP_NUM_THREADS=$NSLOTS
# echo "Got $NSLOTS processors."
. /home/grace/.gmx

#set -u
set -e
set -x

NPME=-1
sysName=sys${SGE_TASK_ID}
NRESUBMITS=2
num=$NUM
base_dir=$PWD

function run {
	MAXH=$1
	cpt_file=sys${SGE_TASK_ID}_prod.cpt
	cd $base_dir/${SGE_TASK_ID}
	echo 'DEBUG: starting mdrun for $MAXH'

	mpirun mdrun -s ${sysName}_prod -deffnm ${sysName}_prod -maxh $MAXH -cpt 720 -nosum -dlb auto -npme $NPME -cpo ${sysName}_prod -cpi $cpt_file
}

run 46

if [ "$num" -lt "$NRESUBMITS" ]; then
	num=$((NUM+1))
	echo "resubmitting - sequence $num for replica $SGE_TASK_ID"
	qsub -v NUM=$num,SGE_TASK_ID=${SGE_TASK_ID} -N ${JOB_NAME}_${SGE_TASK_ID}_${num} /home/grace/labwork/abeta_analysis/scripts/rerun_colosse.sh
fi
