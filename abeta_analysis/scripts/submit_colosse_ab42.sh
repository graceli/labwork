#!/bin/bash
#$ -N abeta42_glca 
#$ -P uix-840-ac
#$ -A uix-840-ac
#$ -l h_rt=48:00:00
#$ -pe default 48 
#$ -q med
#$ -S /bin/bash
#$ -cwd
#$ -notify

module load compilers/intel/11.1.059
module load mpi/openmpi/1.3.4_intel
export OMP_NUM_THREADS=$NSLOTS
# echo "Got $NSLOTS processors."
. /home/grace/.gmx

#set -e
set -u
set -x #trace the bash script

#my variables
num=$NUM #sequence number of the run
NRESUBMITS=2
echo "num is $num"
echo "NRESUBMITS is $NRESUBMITS"

#SHM="/dev/shm/grace"
base_dir=$PWD
echo "starting run in $PWD"

function clean {
        cp -rp $SHM/{*.xtc,*.tpr,*.log,*.cpt,*.edr} $base_dir/sys${SGE_TASK_ID}
	cd $base_dir
        rm -rf $SHM
}

#-notify sends SIGUSR2 
#trap "clean; exit 0" KILL INT SIGINT TERM EXIT SIGUSR2

cd sys$SGE_TASK_ID

cptfile=""
if [ "$num" == "0" ]; then
	grompp -f ../em.mdp -c sys${SGE_TASK_ID}_start.gro -p sys${SGE_TASK_ID}.top -o em.tpr
	mdrun -s em.tpr -deffnm em
	grompp -f ../production.mdp -c em.gro -p sys${SGE_TASK_ID}.top -o prod.tpr
	#cp -p *.tpr $SHM
else
	#cp -p *.cpt *.tpr $SHM
	cptfile="$base_dir/sys${SGE_TASK_ID}/sys${SGE_TASK_ID}_prod.cpt"
fi

#cd $SHM
mpirun /software/apps/gromacs-4.0.7/bin/mdrun -deffnm sys${SGE_TASK_ID}_prod -s $base_dir/sys${SGE_TASK_ID}/prod.tpr -nosum -dlb yes -npme -1 -cpi $cptfile -cpt 60 -maxh 47

# resubmit
#if [ "$num" -lt "$NRESUBMITS" ]; then
#    num=$(($num+1))
#    echo "resubmitting - sequence $num for replica $SGE_TASK_ID"
#    ssh colosse1 "cd $base_dir; qsub -v NUM=$num,SGE_TASK_ID=${SGE_TASK_ID} -N ${JOB_NAME}_${SGE_TASK_ID}_${num} ./run.sh"
#fi

