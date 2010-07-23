#!/bin/bash
#$ -v LD_LIBRARY_PATH=/usr/lib64/tls:/hpf/tools/n1ge/lib/lx24-amd64:/tools/gcc/4.1.1/lib64:/tools/openmpi/1.2/lib:/tools/local/lib
#$ -v PATH=/hpf/tools/n1ge/bin/lx24-amd64:/usr/local/bin:/usr/bin:/usr/X11R6/bin:/bin:/usr/games:/opt/gnome/bin:/tools/local/bin:/home/cneale/exe:/projects/pomes/cneale/exe/gromacs-3.3.1/exec/fftw-3.1.2/bin:.:/tools/openmpi/1.2/bin

### for openmpi LD_LIBRARY_PATH includes /tools/openmpi/1.2b3/lib
### for lam-mpi LD_LIBRARY_PATH includes /tools/lam/lam-7.1.2/lib
### for openmpi PATH includes /tools/openmpi/1.2b3/bin
### for lam-mpi PATH includes /tools/lam/lam-7.1.2/bin

MYMOL=SAME_AS_MYMOL_IN_DRIVER
ED=/projects/pomes/cneale/exe/gromacs-3.3.1/exec/fftw-3.1.2/bin
OMPI=/tools/openmpi/1.2/bin
LAM=/tools/lam/lam-7.1.2/bin
MD=SAME_AS_MD_IN_DRIVER
cd ${MD}

in="_in"
out="_out"
mynp=_mynp

# use scratch directory to avoid NFS delay file corruption problems
cd ${TMPDIR}
cp ${MD}/${MYMOL}${out}.tpr $TMPDIR


#echo "Starting mdrun..."
#Production dynamics
### for openmpi ${OMPI}/mpirun ${ED}/mdrun_openmpi
### for lam-mpi ${LAM}/mpirun C ${ED}/mdrun_mpi

${OMPI}/mpirun ${ED}/mdrun_openmpi_v1.2 -np ${mynp} -s ${MYMOL}${out}.tpr -deffnm ${MYMOL}${out} > output.${MYMOL}_mdrun${out} 2> output.${MYMOL}_mdrun${out}_e
#echo "mdrun finished"

#go to the run directory on the submit node
cd ${MD}

#move the finished simulation files in TMPDIR to MD
mv ${TMPDIR}/* .
for((;;));do
  anyleft=`ls -l ${TMPDIR} | head -1 | awk '{print $2}'`
  if((anyleft==0)); then break; fi
  sleep 5
done

#when moving is done, create a dummy "finished" file to signal run completion
echo ${out} > ./finished

###EOF
