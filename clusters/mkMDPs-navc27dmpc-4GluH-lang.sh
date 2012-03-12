#!/bin/bash

#+++++++++++++++++++++++++++
# for id0 in  1 ; do 
  for((id0=2;id0<=10;id0++)); do
#+++++++++++++++++++++++++++

  base=lang4h-n
  base1=lang4h-n
  r=`echo ${RANDOM}` 
  f="$base""$id0"
  f1="$base1""$id0"

#+++++++++++++++++++++++++++
# z=`awk 'NR=='$id0' {print $5}' IONS.gro`
#+++++++++++++++++++++++++++

  cat langevin_dmpc.mdp_H               > $f.mdp 
  echo "gen_seed          = $r  " >> $f.mdp
#+++++++++++++++++++++++++++

# cat sd.mdp_T                    >> $f.mdp 
# d=dir
# folder="$d""$id0"
# mkdir ./${folder}
# cat $f.mdp > ./${folder}/$f.mdp

#++++++++++++++++++++++++++++ job ++++++++++++++++++++++++++++++++++
  echo \#\!/bin/bash    > $f1.job 
  echo \#PBS -q qwork@mp2 -l walltime=48:00:00 -l nodes=1:ppn=24 >>  $f1.job
  echo \#PBS -N nav4hds$id0      >>  $f1.job
  echo >>  $f1.job
  echo \# to submit this type bqsub eth.sh >>  $f1.job
  echo >>  $f1.job
  echo if \[ \"\$PBS_ENVIRONMENT\" \!= \"PBS_INTERACTIVE\" \]\; then >>  $f1.job
  echo \  if \[ -n \"\$PBS_O_WORKDIR\" \]\; then  >>  $f1.job
  echo \  \  cd \$PBS_O_WORKDIR >>  $f1.job
  echo \  fi  >>  $f1.job
  echo fi >>  $f1.job
  echo >>  $f1.job
  echo . /home/chakraba/.gmx4.5.4_cn >>  $f1.job
  echo sleep 1 >>  $f1.job
  echo >>  $f1.job
#### grompp part ########### only in .job NOT in .rest-job
#### grompp USING gmx454 ###
  echo grompp -f $f.mdp -c t0_4gluh_nav_dmpc_salt.gro  -o $f.tpr -p nilu_4GluH_dmpc_salt.top -po mdout.$f.mdp  >>  $f1.job
  echo sleep 3 >>  $f1.job
  echo module purge >>  $f1.job
  echo module load intel64/12.0.5.220 openmpi_intel64/1.4.3_ofed >>  $f1.job
  echo sleep 5 >>  $f1.job
  echo >>  $f1.job
#### now doing Langevin Dynamics ####
  echo mpirun -np 24 -machinefile \$PBS_NODEFILE2 mdrun_mpi -deffnm $f -s $f.tpr -dlb yes -npme -1 -cpt 60 -maxh 47.8 -cpi $f.cpt -noappend >>  $f1.job
  echo >>  $f1.job
  echo \#\#EOF >>  $f1.job
  chmod u+x $f1.job

#++++++++++++++++++++ rest-job ++++++++++++++++++++++++++++++++++++++++++
  echo \#\!/bin/bash    > $f1.rest-job 
  echo \#PBS -q qwork@mp2 -l walltime=48:00:00 -l nodes=1:ppn=24 >>  $f1.rest-job
  echo \#PBS -N nav4hds$id0      >>  $f1.rest-job
  echo >>  $f1.rest-job
  echo \# to submit this type bqsub eth.sh >>  $f1.rest-job
  echo >>  $f1.rest-job
  echo if \[ \"\$PBS_ENVIRONMENT\" \!= \"PBS_INTERACTIVE\" \]\; then >>  $f1.rest-job
  echo \  if \[ -n \"\$PBS_O_WORKDIR\" \]\; then  >>  $f1.rest-job
  echo \  \  cd \$PBS_O_WORKDIR >>  $f1.rest-job
  echo \  fi  >>  $f1.rest-job
  echo fi >>  $f1.rest-job
  echo >>  $f1.rest-job
  echo . /home/chakraba/.gmx4.5.4_cn >>  $f1.rest-job
  echo >>  $f1.rest-job
  echo sleep 1 >>  $f1.rest-job
  echo module purge >>  $f1.rest-job
  echo module load intel64/12.0.5.220 openmpi_intel64/1.4.3_ofed >>  $f1.rest-job
  echo sleep 5 >>  $f1.job
  echo >>  $f1.rest-job
#### grompp part ########### only in .job NOT in .rest-job ####
#### grompp USING gmx454 ###
#### now doing Langevin Dynamics ####
  echo mpirun -np 24 -machinefile \$PBS_NODEFILE2 mdrun_mpi -deffnm $f -s $f.tpr -dlb yes -npme -1 -cpt 60 -maxh 47.8 -cpi $f.cpt -noappend >>  $f1.rest-job
  echo >>  $f1.rest-job
  echo \#\#EOF >>  $f1.rest-job
  chmod u+x $f1.rest-job
#++++++++++++++++++++ rest-job ++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++
##### don't submit *.rest-job here do that in a restart script ####
#+++++++++++++++++++++++++++

  sleep 2
  bqsub $f1.job

#+++++++++++++++++++++++++++
#

  done
########### EOF - nilu ccccccccccccccccccccccc
