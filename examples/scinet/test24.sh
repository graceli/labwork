#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=24:00:00
#PBS -N gnnq_test
#test per-node efficiency
#mdrun -s gnnq_8_popc_2_r_md1.tpr -deffnm oneproc/gnnq_oneproc
cd /scratch/nikolia/test
mkdir fourproc_x2
echo $PBS_NODEFILE
head -n 4 ${PBS_NODEFILE}
/scinet/gpc/mpi/openmpi/1.3.2-intel-v11.0-ofed/bin/mpirun -np 4 -machinefile $PBS_NODEFILE /scratch/cneale/GPC/exe/intel/gromacs-4.0.5/exec/bin/mdrun_openmpi -deffnm fourproc_x2/gnnq_test_a -s gnnq_8_popc_2_r_md1.tpr -nosum -dlb yes -npme -1 &> output.job1 &
/scinet/gpc/mpi/openmpi/1.3.2-intel-v11.0-ofed/bin/mpirun -np 4 -machinefile $PBS_NODEFILE /scratch/cneale/GPC/exe/intel/gromacs-4.0.5/exec/bin/mdrun_openmpi -deffnm fourproc_x2/gnnq_test_b -s gnnq_8_popc_2_r_md1.tpr -nosum -dlb yes -npme -1 &> output.job2 &
wait
