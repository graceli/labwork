#!/bin/sh 

set -x
set -e

cd $PBS_O_WORKDIR

for xtc in `ls *.xtc`; do
    bname=`basename $xtc .xtc`

    # create a smaller trajectory to make the volume maps faster
    trjconv -f $bname -o ${bname}_smaller.xtc -skip 5
done
