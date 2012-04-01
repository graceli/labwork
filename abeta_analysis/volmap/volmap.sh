#!/bin/sh

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
    trjconv -f /dev/shm/grace/$xtc_file -s $TPR -fit rot+trans -o ${xtc_file}_fit2.xtc
    rm /dev/shm/grace/*
}

mkdir /dev/shm/grace

cd $PBS_O_WORKDIR

is_cer=0

for ratio in 15 64; do 
	for iso in glycerol chiro scyllo; do
        DATA="../$ratio/${iso}_nonsolvent"
		xtc_file="${iso}_${ratio}_volmap_all"
		TPR="../${iso}_${ratio}_nosol.tpr"
		
		concatenate_and_fit $DATA $TPR $xtc_file
		
		if [ "$iso" == "glycerol" ]; then
			is_cer="1"
		fi

        # Spawn a new node for calculating the volume map
		ssh gpc01 "cd $PBS_O_WORKDIR; qsub -v XTC=${xtc_file}_fit.xtc,GRO=${iso}_${ratio}_nosol.gro,is_cer=$is_cer volmap_run.sh"
	done
done
