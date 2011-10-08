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

#-notify sends SIGUSR2 
#trap "clean; exit 0" KILL INT SIGINT TERM EXIT SIGUSR2

if [ ! -e "sys${SGE_TASK_ID}"]; then
	mkdir sys$SGE_TASK_ID
fi
cd sys$SGE_TASK_ID

cptfile=""
if [ "$num" == "0" ]; then
	if [[ ! -e "sys${SGE_TASK_ID}_prod.cpt" || ! -s "sys${SGE_TASK_ID}_prod.cpt" ]]; then	
		grompp -f $base_dir/equil_nvt.mdp -c em.gro -p $base_dir/abeta42.top -o sys${SGE_TASK_ID}_equil_nvt.tpr 
		mpirun mdrun -deffnm sys${SGE_TASK_ID}_equil_nvt -s sys${SGE_TASK_ID}_equil_nvt.tpr -nosum -dlb yes -npme -1 -maxh 47
		
		grompp -f $base_dir/equil_npt.mdp -c sys${SGE_TASK_ID}_equil_nvt -p $base_dir/abeta42.top -o sys${SGE_TASK_ID}_equil_npt
		mpirun mdrun -deffnm sys${SGE_TASK_ID}_equil_npt -s sys${SGE_TASK_ID}_equil_npt -nosum -dlb yes -npme -1 -maxh 47

		grompp -f $base_dir/abeta42_prod.mdp -c sys${SGE_TASK_ID}_equil_npt -p $base_dir/abeta42.top -o sys${SGE_TASK_ID}_prod
	fi
else
	cptfile="$base_dir/${SGE_TASK_ID}/sys${SGE_TASK_ID}_prod.cpt"
fi

mpirun mdrun -deffnm sys${SGE_TASK_ID}_prod -s $base_dir/sys${SGE_TASK_ID}/prod.tpr -nosum -dlb yes -npme -1 -cpi $cptfile -cpt 60 -maxh 47

# resubmit
# if [ "$num" -lt "$NRESUBMITS" ]; then
#    num=$(($num+1))
#    echo "resubmitting - sequence $num for replica $SGE_TASK_ID"
#    ssh colosse1 "cd $base_dir; qsub -v NUM=$num,SGE_TASK_ID=${SGE_TASK_ID} -N ${JOB_NAME}_${SGE_TASK_ID}_${num} ./run.sh"
# fi

