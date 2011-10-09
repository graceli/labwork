#!/bin/bash
#PBS -l nodes=2:ib:ppn=8,walltime=24:00:00
#PBS -N gnnq_test
#test per-node efficiency
#mdrun -s gnnq_8_popc_2_r_md1.tpr -deffnm oneproc/gnnq_oneproc
cd /scratch/nikolia/test
mkdir 16proc
/scinet/gpc/mpi/openmpi/1.3.2-intel-v11.0-ofed/bin/mpirun -np 16 -machinefile $PBS_NODEFILE /scratch/cneale/GPC/exe/intel/gromacs-4.0.5/exec/bin/mdrun_openmpi -deffnm 16proc/gnnq_test -s gnnq_8_popc_2_r_md1.tpr -nosum -dlb yes -npme -1

