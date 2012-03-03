#!/bin/bash
#PBS -l nodes=7:ppn=8,walltime=48:00:00 
#PBS -N abeta42_glca

#export OMP_NUM_THREADS=$NSLOTS
. /home/grace/.gmx

#set -e
set -u
set -x #trace the bash script

#my variables
num=$NUM #sequence number of the run
NRESUBMITS=2
EXE=/global/software/gromacs/gromacs407_intel111/bin/mdrun_d_mpi

echo "num is $num"
echo "NRESUBMITS is $NRESUBMITS"

cd $PBS_O_WORKDIR

base_dir=$PWD

echo "starting run in $PWD"

#-notify sends SIGUSR2 
#trap "clean; exit 0" KILL INT SIGINT TERM EXIT SIGUSR2

if [ ! -e "sys${PBS_ARRAYID}" ]; then
	mkdir sys$PBS_ARRAYID
	cp em.gro sys$PBS_ARRAYID
fi
cd sys$PBS_ARRAYID

cptfile=""
if [ "$num" == "0" ]; then
	if [[ ! -e "sys${PBS_ARRAYID}_prod.cpt" || ! -s "sys${PBS_ARRAYID}_prod.cpt" ]]; then	
		grompp -f $base_dir/equil_nvt.mdp -c em.gro -p $base_dir/abeta42_glucose.top -o sys${PBS_ARRAYID}_equil_nvt.tpr 
		mpirun $EXE -deffnm sys${PBS_ARRAYID}_equil_nvt -s sys${PBS_ARRAYID}_equil_nvt.tpr -nosum -dlb yes -npme -1 -maxh 47
		
		grompp -f $base_dir/equil_npt.mdp -c sys${PBS_ARRAYID}_equil_nvt -p $base_dir/abeta42_glucose.top -o sys${PBS_ARRAYID}_equil_npt
		mpirun $EXE -deffnm sys${PBS_ARRAYID}_equil_npt -s sys${PBS_ARRAYID}_equil_npt -nosum -dlb yes -npme -1 -maxh 47

		grompp -f $base_dir/abeta42_prod.mdp -c sys${PBS_ARRAYID}_equil_npt -p $base_dir/abeta42_glucose.top -o sys${PBS_ARRAYID}_prod
	fi
else
	cptfile="$base_dir/sys${PBS_ARRAYID}/sys${PBS_ARRAYID}_prod.cpt"
fi

if [ ! -e "sys${PBS_ARRAYID}_prod.tpr" ]; then
	grompp -f $base_dir/abeta42_prod.mdp -c sys${PBS_ARRAYID}_equil_npt -p $base_dir/abeta42_glucose.top -o sys${PBS_ARRAYID}_prod
fi

mpirun $EXE -deffnm sys${PBS_ARRAYID}_prod -s $base_dir/sys${PBS_ARRAYID}/sys${PBS_ARRAYID}_prod.tpr -nosum -dlb yes -npme -1 -cpi $cptfile -cpt 60 -maxh 47

# resubmit
#if [ "$num" -lt "$NRESUBMITS" ]; then
#    num=$(($num+1))
#    echo "resubmitting - sequence $num for replica $PBS_ARRAYID"
#   qsub -v NUM=$num,PBS_ARRAYID=${PBS_ARRAYID} -N ${JOB_NAME}_${PBS_ARRAYID}_${num} /home/grace/labwork/abeta_analysis/scripts/submit_colosse_ab42.sh
#fi

