#!/bin/sh
#PBS -l nodes=1:compute-eth:ppn=8,walltime=02:00:00,os=centos53computeA
#PBS -N analysis

#set -u
#set -e
set -x

proteingrp=1
insgrp=12
chain1=18
chain2=19
chain3=20
chain4=21
chain5=22
xvgr="-noxvgr"
start=1
end=20

#extract the number of nonpolar contacts between inositol and residues
# function make_indices {
# 	echo -e "'SideChain'&aC*&!rACE\nsplitch17\nq" | make_ndx -f common/em.tpr -o common/ab_nonpolar.ndx 
# 	RETURN_CODE=$?
# }

# function nonpolar {
# 	for s in `seq $start $end`; do
# 		DATA=../sys${s}
# 		#mkdir sys${s}_o
# 		echo -e "'SideChain'&aC*&!rACE\nsplitch17\nq" | make_ndx -f $DATA/em.tpr -o $DATA/ab_nonpolar.ndx 
# 		for xtc in $DATA/sys${s}.xtc $DATA/sys${s}.part0002.xtc $DATA/sys${s}.part0003.xtc; do 
# 			echo "analyzing $xtc ..."
# 			for i in $chain1 $chain2 $chain3 $chain4 $chain5; do
# 				echo $i $insgrp | g_inositol_residue_nonpolar_v2 -f $xtc -s $DATA/em.tpr -n $DATA/ab_nonpolar.ndx -per_residue_contacts chain${i}_residue_np_contact.dat -per_inositol_contacts chain${i}_inositol_np_contact.dat
# 				 mv table.dat chain${i}_table.dat
# 			
# 				trap "exit 1" SIGTERM TERM KILL SIGINT
# 			done
# 
# 			### evoke python script ###
# 			python abeta_analysis.py sys${s}.h5
# 			###########################
# 			rm *.dat
# 		done
# 	done
# }

# function dssp {
#         for ratio in 15 64; do
#                 for iso in scyllo chiro; do
#                         trj="ab_${ratio}_${iso}.xtc_c_fit.xtc"
#                         echo 1 | do_dssp -f $ratio/$trj -s $ratio/nosol.tpr -o ab_${ratio}_${iso}_ss -sc ab_${ratio}_${iso}_sc -dt 10
#                         trap "exit 1" TERM INT KILL
#                 done
#         done
# }

trap "clean; exit" TERM INT SIGINT

function clean {
	cd /dev/shm
	tar cvfz "analysis_${1}.tgz" analysis
	cp analysis*.tgz $base_dir
	rm -rf analysis analysis.tgz
}

# calculate the rmsd of the protein using the nmr structure as a reference
function rmsd {
    iso=$1
    ratio=$2
	output_dir=$3/rmsd
	mkdir -p $output_dir
	
	seq 1 10 | parallel -j 8 "echo 1 1 | g_rms -f $DATA/ab_${iso}_${ratio}_{}_nosol_whole.xtc_c_fit.xtc -s nmr_protein.tpr -o $output_dir/ab_${iso}_${ratio}_{}_nosol_whole_rmsd_protein.xvg"
	seq 1 10 | parallel -j 8 "echo 4 4 | g_rms -f $DATA/ab_${iso}_${ratio}_{}_nosol_whole.xtc_c_fit.xtc -s nmr_protein.tpr -o $output_dir/ab_${iso}_${ratio}_{}_nosol_whole_rmsd_backbone.xvg"
	
	clean "${iso}_${ratio}_rmsd"
}

# calculate the rmsf
function rmsf_calpha {
	iso=$1
    ratio=$2
	output_dir=$3/rmsd
	mkdir -p $output_dir
	
    calpha=3
    # backbone fitting for specific parts of the peptide    
	seq 1 10 | parallel -j 8 "echo $calpha | g_rmsf -f ab_${iso}_${ratio}_${i}_nosol_whole.xtc_c_fit.xtc -s ${ratio}_nosol.tpr -o $output_dir/${iso}_${ratio}_rmsf.xvg -fit -res -ox $output_dir/ab_${iso}_${ratio}_${i}_nosol_whole.xtc_c_fit -noxvgr"
	clean "${iso}_${ratio}_rmsf"
}

cd $PBS_O_WORKDIR
base_dir=`pwd`
DATA="/scratch/grace/inositol/abeta42/2/xtc"
SHM="/dev/shm/analysis"

#Externally bound variables
#ISO=scyllo
#RATIO=15
#ANALYSIS=rmsd
echo $ISO
echo $RATIO
echo $ANALYSIS

${ANALYSIS} $ISO $RATIO $SHM


