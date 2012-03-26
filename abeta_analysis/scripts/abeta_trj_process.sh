#!/bin/bash

#PBS -l nodes=1:ppn=8,walltime=2:00:00,os=centos53computeA
#PBS -N abeta_process

set -u
#set -e 
set -x 

trap "clean_exit; exit 0" TERM INT KILL SIGINT SIGTERM EXIT

SHM=/dev/shm/grace
mkdir -p $SHM

function clean_exit {
	cd $PBS_O_WORKDIR
	output_persistent="xtc" #$RATIO/$ISO/$REPLICA
	mkdir -p $output_persistent
		
	cd $SHM
	cp -rp . $PBS_O_WORKDIR/$output_persistent
	rm -rf *
}

protein_ins=1
protein=0
center=2
function center {
        trj=$1
		tpr=$2
		ndx=$3
		input_dir=$4	
        #centering
        echo $center $protein_ins | trjconv -f $input_dir/whole/$trj -s $tpr -center -pbc res -o $output_dir/${trj}_c.xtc -n $ndx
        #fitting
        echo $protein $protein_ins | trjconv -f $output_dir/${trj}_c.xtc -s $tpr -fit rot+trans -o $output_dir/${trj}_c_fit.xtc -n $ndx
		rm  $output_dir/${trj}_c.xtc
}


cd $PBS_O_WORKDIR

input_dir="$RATIO/$ISO/array/"
output_dir=$SHM

#trjcat
#echo 1 | trjcat -f $input_dir/$REPLICA/*prod*.xtc -o $output_dir/ab_${ISO}_${RATIO}_${REPLICA}_nosol.xtc -n ${ISO}_${RATIO}_nosol.ndx

#trjconv
#echo 1 | trjconv -f $output_dir/ab_${ISO}_${RATIO}_${REPLICA}_nosol.xtc -s ${ISO}_${RATIO}_nosol.tpr -pbc whole -o $output_dir/ab_${ISO}_${RATIO}_${REPLICA}_nosol_whole.xtc -n ${ISO}_${RATIO}_nosol.ndx

#rm $output_dir/ab_${ISO}_${RATIO}_${REPLICA}_nosol.xtc

center ab_${ISO}_${RATIO}_${REPLICA}_nosol_whole.xtc  ${ISO}_${RATIO}_nosol.tpr  ${ISO}_${RATIO}_nosol.ndx xtc
#center ab_15_scyllo.xtc
