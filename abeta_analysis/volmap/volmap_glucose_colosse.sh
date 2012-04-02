#!/bin/sh
#PBS -l nodes=1:ppn=8,walltime=8:00:00
#PBS -N volmap_full

set -x

trap 'clean "${iso}_${r}"; exit $?' TERM KILL EXIT SIGINT

function clean {
    cd /dev/shm/grace
    tar cvfz $PBS_O_WORKDIR/${1}_volmap_full_dt10.tar.gz *
    rm -rf /dev/shm/grace
}

mkdir /dev/shm/grace

cd $PBS_O_WORKDIR

XTC=${iso}_${r}_volmap_all_fit2.xtc
GRO=${iso}_${r}_nosol.gro

if [ "$iso" == "glycerol" ]; then
    vmd -dispdev text -e volmap_cer.tcl -args $XTC $GRO > /dev/shm/grace/${XTC}_volmap_full.log 2>&1 &
else
    vmd -dispdev text -e volmap.tcl -args $XTC $GRO > /dev/shm/grace/${XTC}_volmap_full.log 2>&1 &
fi

wait
#!/bin/sh
#PBS -l nodes=1:ppn=8,walltime=8:00:00
#PBS -N trjcat_for_volmap

# Batch compute "volume maps" of solutes chiro scyllo and glycerol around Abeta42 using vmd

set -x
set -e

trap 'clean; exit $?' KILL TERM SIGINT

function clean {
    rm -rf /dev/shm/grace
}

function concatenate_and_fit {
    DATA=$1
    TPR=$2
    xtc_file=$3
    
    for i in `seq 1 100`; do echo "c"; done | trjcat -f $DATA/*.xtc -o /dev/shm/grace/$xtc_file -cat -keeplast -settime -dt 10
    # fit on protein group, output nonsolvent group.
    echo "Protein Protein_INS" | trjconv -f /dev/shm/grace/$xtc_file -s $TPR -fit rot+trans -o ${xtc_file}_fit2.xtc -n ${iso}_${ratio}_volmap.ndx
    rm /dev/shm/grace/*
}

mkdir /dev/shm/grace

cd $PBS_O_WORKDIR

is_cer=0

# ratio and iso variables have passed in values via qsub

DATA="../$ratio/${iso}_nonsolvent"
xtc_file="${iso}_${ratio}_volmap_all"
TPR="../${iso}_${ratio}_nosol.tpr"

concatenate_and_fit $DATA $TPR $xtc_file

#if [ "$iso" == "glycerol" ]; then
    #is_cer="1"
#fi

#for ratio in 15 64; do
    #for iso in glycerol chiro scyllo; do
        ## Spawn a new node for calculating the volume map
        #ssh gpc01 "cd $PBS_O_WORKDIR; qsub -v XTC=${xtc_file}_fit.xtc,GRO=${iso}_${ratio}_nosol.gro,is_cer=$is_cer volmap_run.sh"
    #done
#done
