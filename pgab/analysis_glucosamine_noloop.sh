#!/bin/bash
#PBS -q qwork@mp2
#PBS -l walltime=5:00:00 -l nodes=1:ppn=1
#PBS -N pgab

# stricter bash -- quits on error and unset variables

. $HOME/.gmx4.5.4

set -u
set -e
set -x

DATA=/mnt/scratch_mp2/pomes/ligrace1/pgab/pgab_glucosamine/pgab_glucosamine_nonsolvent
function rmsd_no_loop {
    for i in `seq 1 13`; do
        g_rmsdist_mpi -f $DATA/${i}_final.xtc -s $DATA/pgab_glucosamine.tpr -dt 10 -o glucosamine_rmsd_backbone_no_loop_${i}.xvg -n protein_noloop.ndx &
    done
    wait
}

#cd $PBS_O_WORKDIR

rmsd_no_loop
