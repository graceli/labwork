#!/bin/bash
#PBS -q qwork@mp2
#PBS -l walltime=100:00:00 -l nodes=3:ppn=1
#PBS -N abeta 

# stricter bash -- quits on error and unset variables
set -u
set -e
set -x

. $HOME/.gmx4.5.4

#set -e
set -u
set -x #trace the bash script

#my variables

EXE=mdrun_mpi
sysName=system_${PBS_ARRAYID}_final
SHM=/dev/shm
#base_dir=$PBS_O_WORKDIR
NRESUBMITS=0
#num=$NUM
DEBUG=0
MAXH=100
nodes=3
NCORES=$(echo "$nodes*24" | bc)
NPME=24
PARAMS="../params"
EQUIL_NVT_MDP="equil_nvt.mdp"
EQUIL_NPT_MDP="equil_berendsen.mdp"
PROD_MDP="prod_npt.mdp"

num=$NUM #sequence number of the run
NRESUBMITS=2
echo "num is $num"
echo "NRESUBMITS is $NRESUBMITS"

#SHM="/dev/shm/grace"
# Hard coding the base_dir because it was picking up a weird path
base_dir=/mnt/scratch_mp2/pomes/ligrace1/abeta/15/glucose
echo "starting run in /mnt/scratch_mp2/pomes/ligrace1/abeta/15/glucose"

#-notify sends SIGUSR2 
#trap "clean; exit 0" KILL INT SIGINT TERM EXIT SIGUSR2
cd $PBS_O_WORKDIR
cd sys$PBS_ARRAYID

cptfile=""
if [ "$num" == "0" ]; then
	if [[ ! -e "sys${PBS_ARRAYID}_prod.cpt" || ! -s "sys${PBS_ARRAYID}_prod.cpt" ]]; then	
		if [[ ! -e "em.gro" ]]; then
			grompp_mpi -f $base_dir/em.mdp -c ions_counter.gro -p $base_dir/abeta42_glucose_15.top -o ${PBS_ARRAYID}_em.tpr
			mdrun_mpi -deffnm em -s ${PBS_ARRAYID}_em.tpr
		fi

		grompp_mpi -f $base_dir/equil_nvt.mdp -c em.gro -p $base_dir/abeta42_glucose_15.top -o sys${PBS_ARRAYID}_equil_nvt.tpr 
		mpirun mdrun_mpi -deffnm sys${PBS_ARRAYID}_equil_nvt -s sys${PBS_ARRAYID}_equil_nvt.tpr -nosum -dlb yes -npme -1 -maxh 47
		
		grompp_mpi -f $base_dir/equil_npt.mdp -c sys${PBS_ARRAYID}_equil_nvt -p $base_dir/abeta42_glucose_15.top -o sys${PBS_ARRAYID}_equil_npt
		mpirun mdrun_mpi -deffnm sys${PBS_ARRAYID}_equil_npt -s sys${PBS_ARRAYID}_equil_npt -nosum -dlb yes -npme -1 -maxh 47

		grompp_mpi -f $base_dir/abeta42_prod.mdp -c sys${PBS_ARRAYID}_equil_npt -p $base_dir/abeta42_glucose_15.top -o sys${PBS_ARRAYID}_prod
	fi
else
	cptfile="$base_dir/sys${PBS_ARRAYID}/sys${PBS_ARRAYID}_prod.cpt"
fi

if [ ! -e "sys${PBS_ARRAYID}_prod.tpr" ]; then
	grompp_mpi -f $base_dir/abeta42_prod.mdp -c sys${PBS_ARRAYID}_equil_npt -p $base_dir/abeta42_glucose_15.top -o sys${PBS_ARRAYID}_prod
fi

mpirun -np $NCORES $EXE -s $base_dir/sys${PBS_ARRAYID}/sys${PBS_ARRAYID}_prod.tpr -deffnm sys${PBS_ARRAYID}_prod -maxh $MAXH -cpt 720 -nosum -dlb auto -npme $NPME -cpo prod -cpi -noappend

# resubmit
#if [ "$num" -lt "$NRESUBMITS" ]; then
#    num=$(($num+1))
#    echo "resubmitting - sequence $num for replica $PBS_ARRAYID"
#   qsub -v NUM=$num,PBS_ARRAYID=${PBS_ARRAYID} -N ${JOB_NAME}_${PBS_ARRAYID}_${num} /home/grace/labwork/abeta_analysis/scripts/submit_colosse_ab42.sh
#fi

