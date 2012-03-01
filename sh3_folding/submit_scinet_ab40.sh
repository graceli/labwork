#!/bin/bash
#MOAB/Torque submission script for multiple serial jobs on SciNet GPC
#PBS -l nodes=1:ppn=8,walltime=48:00:00,os=centos53computeA
#PBS -N ab40_2fold_1

REPEAT=2
NPME=-1
EXE=/project/pomes/cneale/GPC/exe/intel/gromacs-4.0.5/exec/bin/mdrun_openmpi-1.4.1
sysName=sys1
top=abeta40_water.top
cd $PBS_O_WORKDIR

#grompp -f em.mdp -c ions_counter.gro -p abeta42.top -o em.tpr
mdrun -s em.tpr -deffnm em

grompp -f equil_nvt.mdp -c em.gro -p $top -o ${sysName}_equil_nvt.tpr 
/scinet/gpc/mpi/openmpi/1.4.1-intel-v11.0-ofed/bin/mpirun -np $(wc -l $PBS_NODEFILE | gawk '{print $1}') -machinefile $PBS_NODEFILE $EXE -s ${sysName}_equil_nvt.tpr -deffnm ${sysName}_equil_nvt -npme $NPME

grompp -f equil_npt.mdp -c ${sysName}_equil_nvt -p $top -o ${sysName}_equil_npt
/scinet/gpc/mpi/openmpi/1.4.1-intel-v11.0-ofed/bin/mpirun -np $(wc -l $PBS_NODEFILE | gawk '{print $1}') -machinefile $PBS_NODEFILE $EXE -s ${sysName}_equil_npt.tpr -deffnm ${sysName}_equil_npt -npme $NPME 

grompp -f abeta42_prod.mdp -c em.gro -p $top -o ${sysName}_prod.tpr
/scinet/gpc/mpi/openmpi/1.4.1-intel-v11.0-ofed/bin/mpirun -np $(wc -l $PBS_NODEFILE | gawk '{print $1}') -machinefile $PBS_NODEFILE $EXE  -s ${sysName}_prod -deffnm ${sysName}_prod -maxh 47 -cpt 720 -nosum -dlb yes -npme $NPME  -cpo ${sysName}_prod -cpi 


