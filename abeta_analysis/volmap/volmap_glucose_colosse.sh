#!/bin/sh
# $ -N ab_volmap

set -e
set -x

function clean {
    cd /dev/shm/grace
    tar cvfz volmap_glucose.tar.gz *
    rm -rf /dev/shm/grace
}

function concatenate_and_fit {
    DATA=$1
    TPR=$2
    xtc_file=$3
    
    for i in `seq 1 100`; do echo "c"; done | trjcat -f $DATA/*.xtc -o /dev/shm/grace/$xtc_file -cat -keeplast -settime -dt 10
    
    # fit on protein group, output nonsolvent group.
    echo "Protein Protein_GLCA" | trjconv -f /dev/shm/grace/$xtc_file -s $TPR -fit rot+trans -o ${xtc_file}_fit2.xtc -n protein_glca.ndx
    
    rm /dev/shm/grace/*
}

mkdir /dev/shm/grace

# ratio and iso variables have passed in values via qsub
# Batch compute "volume maps" of glucose around Abeta

DATA="../../${iso}_Protein_GLCA"
xtc_file="${iso}_volmap_all"
TPR="protein_glca.tpr"
GRO="protein_glca.gro"

concatenate_and_fit $DATA $TPR $xtc_file

vmd -dispdev text -e /home/grace/labwork/abeta_analysis/volmap/volmap_glca.tcl -args $xtc_file $GRO > /dev/shm/grace/${xtc_file}.log 2>&1 &

