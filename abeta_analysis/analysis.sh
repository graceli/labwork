#!/bin/sh
#PBS -l nodes=1:compute-eth:ppn=8,walltime=05:00:00,os=centos53computeA
#PBS -N analysis

set -u
#set -e
set -x



#extract the number of nonpolar contacts between inositol and residues
# function make_indices {
# 	echo -e "'SideChain'&aC*&!rACE\nsplitch17\nq" | make_ndx -f common/em.tpr -o common/ab_nonpolar.ndx 
# 	RETURN_CODE=$?
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

proteingrp=1
insgrp=12
chain1=17
chain2=18
chain3=19
chain4=20
chain5=21
xvgr="-noxvgr"
# start=1
# end=20
function nonpolar {
	iso=$1
    ratio=$2
	output_dir=$3/rmsd
	mkdir -p $output_dir
	
	# echo -e "'SideChain'&aC*&!rACE\nsplitch17\nq" | make_ndx -f ${ratio}_nosol.tpr -o ab_nonpolar.ndx
	make_ndx -f ${ratio}_nosol.tpr -o ab_nonpolar.ndx << EOF
	'SideChain'&aC*&!rACE
	splitch16
	q
	EOF
	
	for s in `seq 1 10`; do
		xtc="ab_${iso}_${ratio}_${s}_nosol_whole.xtc_c_fit"
		if [ -e "$DATA/$xtc" ]; then
			seq $chain1 $chain5 | parallel -j 5 "echo $i $insgrp | g_inositol_residue_nonpolar_v2 -f $DATA/$xtc -s ${ratio}_nosol.tpr -n ab_nonpolar.ndx -per_residue_contacts $output_dir/chain${i}_residue_np_contact.dat -per_inositol_contacts $output_dir/chain${i}_inositol_np_contact.dat"
		fi
	done
}

# calculate the rmsd of the protein using the nmr structure as a reference
function rmsd {
    iso=$1
    ratio=$2
	output_dir=$3/rmsd
	mkdir -p $output_dir
	
	seq 1 10 | parallel -j 8 "echo 1 1 | g_rms -f $DATA/ab_${iso}_${ratio}_{}_nosol_whole.xtc_c_fit.xtc -s nmr_protein.tpr -o $output_dir/ab_${iso}_${ratio}_{}_nosol_whole_rmsd_protein.xvg -noxvgr"
	seq 1 10 | parallel -j 8 "echo 4 4 | g_rms -f $DATA/ab_${iso}_${ratio}_{}_nosol_whole.xtc_c_fit.xtc -s nmr_protein.tpr -o $output_dir/ab_${iso}_${ratio}_{}_nosol_whole_rmsd_backbone.xvg -noxvgr"
	
	clean "${iso}_${ratio}_rmsd"
}

# calculate the rmsf
function rmsf_calpha {
	iso=$1
    ratio=$2
	output_dir=$3/rmsf
	mkdir -p $output_dir
    calpha=3

    # backbone fitting for specific parts of the peptide    
	seq 1 10 | parallel -j 8 "echo $calpha | g_rmsf -f $DATA/ab_${iso}_${ratio}_{}_nosol_whole.xtc_c_fit.xtc -s ${ratio}_nosol.tpr -o $output_dir/${iso}_${ratio}_{}_rmsf.xvg -fit -res -ox $output_dir/ab_${iso}_${ratio}_{}_nosol_whole.xtc_c_fit -noxvgr $TEST" 
	clean "${iso}_${ratio}_rmsf"
}

cd $PBS_O_WORKDIR
base_dir=`pwd`
DATA="/scratch/grace/inositol/abeta42/2/xtc"
SHM="/dev/shm/analysis"

#Externally bound variables
#ISO=water
#RATIO=64
#ANALYSIS=rmsf_calpha
#echo $ISO
#echo $RATIO
#echo $ANALYSIS

TEST="-b 0"
echo `pwd`
${ANALYSIS} $ISO $RATIO $SHM


