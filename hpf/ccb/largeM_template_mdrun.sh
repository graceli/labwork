#!/bin/bash
#$ -v LD_LIBRARY_PATH=/usr/lib64/tls:/hpf/tools/n1ge/lib/lx24-amd64:/tools/gcc/4.1.1/lib64:/tools/openmpi/1.2.1/lib:/tools/local/lib
#$ -v PATH=/hpf/tools/n1ge/bin/lx24-amd64:/usr/local/bin:/usr/bin:/usr/X11R6/bin:/bin:/usr/games:/opt/gnome/bin:/tools/local/bin:/home/cneale/exe:/projects/pomes/cneale/exe/gromacs-3.3.1/exec/fftw-3.1.2/bin:.:/tools/openmpi/1.2.1/bin

### for lam-mpi LD_LIBRARY_PATH includes /tools/lam/lam-7.1.2/lib
### for openmpi_v1.2.1 LD_LIBRARY_PATH includes tools/openmpi/1.2.1/lib and /tools/local/lib
### for lam-mpi PATH includes /tools/lam/lam-7.1.2/bin
### for openmpi_v1.2.1 PATH includes /tools/openmpi/1.2.1/bin

MYMOL=sys225
ED=/projects/pomes/cneale/exe/gromacs-3.3.1/exec/fftw-3.1.2/bin
OMPI=/tools/openmpi/1.2.1/bin
LAM=/tools/lam/lam-7.1.2/bin
MD=/projects/pomes/cneale/micelle/2large
cd ${MD}

in="_in"
out="_out"
mynp=_mynp

echo "Starting mdrun..."

#Production dynamics
### for openmpi ${OMPI}/mpirun ${ED}/mdrun_openmpi_v1.2.1
### for lam-mpi ${LAM}/mpirun C ${ED}/mdrun_mpi

${OMPI}/mpirun ${ED}/mdrun_openmpi_v1.2.1 -np ${mynp} -s ${MYMOL}${out}.tpr -deffnm ${MYMOL}${out} -v > output.${MYMOL}_mdrun${out} 2> output.${MYMOL}_mdrun${out}_e
echo "mdrun finished"

#Deshuffle the gro xtc and trr files. The edr file does not need this.
echo System | ${ED}/trjconv -f ${MYMOL}${out}.xtc -s ${MYMOL}${out}.tpr -n deshuffledesort${MYMOL}${out}.ndx -o ${MYMOL}${out}_deshuffleddesorted.xtc
echo System | ${ED}/trjconv -f ${MYMOL}${out}.trr -s ${MYMOL}${out}.tpr -n deshuffledesort${MYMOL}${out}.ndx -o ${MYMOL}${out}_deshuffleddesorted.trr
echo System | ${ED}/trjconv -f ${MYMOL}${out}.gro -s ${MYMOL}${out}.tpr -n deshuffledesort${MYMOL}${out}.ndx -o ${MYMOL}${out}_deshuffleddesorted.gro


echo ${out} > ./finished

###EOF
