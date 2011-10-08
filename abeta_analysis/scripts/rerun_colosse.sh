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

NPME=-1
sysName=sys${SGE_TASK_ID}
NRESUBMITS=2
num=$NUM
base_dir=$PWD

function run {
	MAXH=$1
	cpt_file=sys${SGE_TASK_ID}_prod.cpt
	cd $base_dir/sys${SGE_TASK_ID}
	echo 'DEBUG: starting mdrun for $MAXH'
	mpirun mdrun -s sys${SGE_TASK_ID}_prod -deffnm sys${SGE_TASK_ID}_prod -cpt 720 -nosum -dlb auto -npme -1 -cpo sys${SGE_TASK_ID}_prod -cpi sys${SGE_TASK_ID}_prod.cpt -maxh 0.03
}

run 46

if [ "$num" -lt "$NRESUBMITS" ]; then
	num=$((NUM+1))
	echo "resubmitting - sequence $num for replica $SGE_TASK_ID"
	qsub -v NUM=$num,SGE_TASK_ID=${SGE_TASK_ID} -N ${JOB_NAME}_${SGE_TASK_ID}_${num} ../run.sh
fi
