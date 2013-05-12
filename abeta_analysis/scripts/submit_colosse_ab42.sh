#!/bin/bash
#PBS -A uix-840-af
#PBS -N ab_glca_15
#PBS -l walltime=48:00:00
#PBS -l nodes=7:ppn=8
#PBS -l flags=restartable
#PBS -l gattr=ckpt
#PBS -o Ab_glca_15-%I.out
#PBS -e Ab_glca_15-%I.err



#path for qsub
export PATH=$PATH:/usr/local/ge6.2u5/bin/lx24-amd64

module load compilers/intel/11.1.059
module load mpi/openmpi/1.4.3_intel
#module load mpi/openmpi/1.3.4_intel
export OMP_NUM_THREADS=$NSLOTS
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

if [ ! -e "sys${MOAB_JOBARRAYINDEX}" ]; then
	mkdir sys$MOAB_JOBARRAYINDEX
	cp em.gro sys$MOAB_JOBARRAYINDEX
fi
cd sys$MOAB_JOBARRAYINDEX

cptfile=""
if [ "$num" == "0" ]; then
	if [[ ! -e "sys${MOAB_JOBARRAYINDEX}_prod.cpt" || ! -s "sys${MOAB_JOBARRAYINDEX}_prod.cpt" ]]; then	
		if [[ ! -e "em.gro" ]]; then
			grompp -f $base_dir/em.mdp -c ions_counter.gro -p $base_dir/abeta42_glucose_15.top -o ${MOAB_JOBARRAYINDEX}_em.tpr
			mdrun -deffnm em -s ${MOAB_JOBARRAYINDEX}_em.tpr
		fi

		grompp -f $base_dir/equil_nvt.mdp -c em.gro -p $base_dir/abeta42_glucose_15.top -o sys${MOAB_JOBARRAYINDEX}_equil_nvt.tpr 
		mpirun mdrun -deffnm sys${MOAB_JOBARRAYINDEX}_equil_nvt -s sys${MOAB_JOBARRAYINDEX}_equil_nvt.tpr -nosum -dlb yes -npme -1 -maxh 47
		
		grompp -f $base_dir/equil_npt.mdp -c sys${MOAB_JOBARRAYINDEX}_equil_nvt -p $base_dir/abeta42_glucose_15.top -o sys${MOAB_JOBARRAYINDEX}_equil_npt
		mpirun mdrun -deffnm sys${MOAB_JOBARRAYINDEX}_equil_npt -s sys${MOAB_JOBARRAYINDEX}_equil_npt -nosum -dlb yes -npme -1 -maxh 47

		grompp -f $base_dir/abeta42_prod.mdp -c sys${MOAB_JOBARRAYINDEX}_equil_npt -p $base_dir/abeta42_glucose_15.top -o sys${MOAB_JOBARRAYINDEX}_prod
	fi
else
	cptfile="$base_dir/sys${MOAB_JOBARRAYINDEX}/sys${MOAB_JOBARRAYINDEX}_prod.cpt"
fi

if [ ! -e "sys${MOAB_JOBARRAYINDEX}_prod.tpr" ]; then
	grompp -f $base_dir/abeta42_prod.mdp -c sys${MOAB_JOBARRAYINDEX}_equil_npt -p $base_dir/abeta42_glucose_15.top -o sys${MOAB_JOBARRAYINDEX}_prod
fi

mpirun /software/apps/gromacs-4.0.7/bin/mdrun -deffnm sys${MOAB_JOBARRAYINDEX}_prod -s $base_dir/sys${MOAB_JOBARRAYINDEX}/sys${MOAB_JOBARRAYINDEX}_prod.tpr -nosum -dlb yes -npme -1 -cpi $cptfile -cpt 60 -maxh 47

# resubmit
#if [ "$num" -lt "$NRESUBMITS" ]; then
#    num=$(($num+1))
#    echo "resubmitting - sequence $num for replica $MOAB_JOBARRAYINDEX"
#   qsub -v NUM=$num,MOAB_JOBARRAYINDEX=${MOAB_JOBARRAYINDEX} -N ${JOB_NAME}_${MOAB_JOBARRAYINDEX}_${num} /home/grace/labwork/abeta_analysis/scripts/submit_colosse_ab42.sh
#fi

