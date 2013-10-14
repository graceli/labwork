#!/bin/bash
#PBS -q qwork@mp2 
#PBS -l walltime=5:00:00 -l nodes=1:ppn=1
#PBS -N pgab

# stricter bash -- quits on error and unset variables

. $HOME/.gmx4.5.4

set -u
set -e
set -x

function rmsd {
	for i in `seq 37 50`; do
		echo Backbone | g_rmsdist_mpi -f ../pgab_nonsolvent/${i}_final.xtc -s ../pgab/pgab.tpr -dt 10 -o rmsd_backbone_${i}.xvg &
	done
	wait
}

function rmsf {
	for i in `seq 37 50`; do
		echo Protein | g_rmsf_mpi -f ../pgab_nonsolvent/${i}_final.xtc -s ../pgab/pgab.tpr -dt 10 -o rmsf_protein_${i}.xvg -od rmsdev_protein_od_${i}.xvg -res &
	done
	wait
}


function rmsf_backbone {
	echo Backbone | g_rmsf_mpi -f ../../volmap/pgab_37-50_all_dt10_fit.xtc -s ../../pgab/pgab.tpr -o pgab_37-50_all_dt10_fit_rmsf_backbone.xvg -od pgab_37-50_all_dt10_fit_rmsdev_backbone_od.xvg -res
}


function rmsd_no_loop {
    for i in `seq 37 50`; do
        g_rmsdist_mpi -f ../pgab_nonsolvent/${i}_final.xtc -s ../pgab/pgab.tpr -dt 10 -o rmsd_backbone_no_loop_${i}.xvg -n protein_no_loop_${i}.ndx &
    done
    wait
}

cd $PBS_O_WORKDIR

rmsd_no_loop

