#!/bin/sh
#PBS -l nodes=1:compute-eth:ppn=8,walltime=02:00:00,os=centos53computeA
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

trap 'echo $?' TERM INT SIGINT EXIT SIGKILL SIGSTOP SIGTERM

function clean {
	cd /dev/shm
	tar cvfz "analysis_${1}.tgz" grace
	cp analysis*.tgz $base_dir
	rm -rf * 
}

function dssp {
	export DSSP=/home/grace/src/dssp_ana/dsspcmbi
	
	for i in `seq 0 10`; do 
		if [ ! -e "$SHM/$i" ]; then
			mkdir -p $SHM/$i
		fi
	done

	seq 0 9 | parallel -j 8 "cd $SHM/{}; echo 1 | do_dssp -f $DATA/ab_${ISO}_${RATIO}_{}_nosol_whole.xtc_c_fit -s $base_dir/${ISO}_${RATIO}_nosol.tpr -o ab_${ISO}_${RATIO}_{}_ss -sc ab_${ISO}_${RATIO}_{}_sc $TEST 2>&1"
	clean "${ISO}_${RATIO}_dssp"
}

chain_start=0
chain_end=3
function chain_hbonds {
	iso=$1
    ratio=$2
	output_dir=$3/hbonds
	mkdir -p $output_dir
		
	for s in `seq 1 10`; do
		xtc="ab_${iso}_${ratio}_${s}_nosol_whole.xtc_c_fit"
		if [ -e "$DATA/${xtc}.xtc" ]; then
			mkdir -p $output_dir/$s
			for ch in `seq $chain_start $chain_end`; do
				let next=ch+1
				echo $ch $next | g_hbond -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n chain.ndx -nonitacc -nomerge -num $output_dir/$s/chain_${ch}_${next}_hbonds $TEST > /dev/null 2>&1 &
			done
			wait
		fi
		# python /home/grace/AnalysisScripts/abeta_analysis/abeta_analysis.py sys${s}.h5
	done
	clean "${iso}_${ratio}_chain_hbonds"
}

res_start=0
res_end=129
INS_grp=130
num=0
function hbonds {
	iso=$1
    ratio=$2
	output_dir=$3/hbonds
	mkdir -p $output_dir
		
	for s in `seq 1 10`; do
		xtc="ab_${iso}_${ratio}_${s}_nosol_whole.xtc_c_fit"
		if [ -e "$DATA/${xtc}.xtc" ]; then
			mkdir -p $output_dir/$s
			seq $res_start $res_end | parallel -j 8 "echo {} $INS_grp | g_hbond -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n g_hbond_${ratio}_${iso}.ndx -nonitacc -nomerge -num $output_dir/$s/{} $xvgr $TEST > /dev/null 2>&1"
		fi
		# python /home/grace/AnalysisScripts/abeta_analysis/abeta_analysis.py sys${s}.h5
	done
	clean "${iso}_${ratio}_hbonds"
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
	output_dir=$3/nonpolar
	mkdir -p $output_dir
	
	 echo -e "'SideChain'&aC*&!rACE\nsplitch16\nq" | make_ndx -f ${iso}_${ratio}_nosol.tpr -o ab_${ratio}_nonpolar.ndx
	
	for s in `seq 1 10`; do
		xtc="ab_${iso}_${ratio}_${s}_nosol_whole.xtc_c_fit"
		if [ -e "$DATA/${xtc}.xtc" ]; then
			seq $chain1 $chain5 | parallel -j 5 "echo {} $insgrp | g_inositol_residue_nonpolar_v2 -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n ab_${ratio}_nonpolar.ndx -per_residue_contacts $output_dir/chain{}_residue_np_contact.dat -per_inositol_contacts $output_dir/chain{}_inositol_np_contact.dat -per_residue_table chain{}_table.dat $TEST"
		fi
	done
	clean "${iso}_${ratio}_nonpolar"
}

# calculate the rmsd of the protein using the nmr structure as a reference
function rmsd {
    iso=$1
    ratio=$2
	output_dir=$3/rmsd
	mkdir -p $output_dir
	
	seq 1 10 | parallel -j 8 "echo 1 1 | g_rms -f $DATA/ab_${iso}_${ratio}_{}_nosol_whole.xtc_c_fit.xtc -s ${iso}_${ratio}_nosol.tpr -o $output_dir/ab_${iso}_${ratio}_{}_nosol_whole_rmsd_protein.xvg -noxvgr $TEST"
	seq 1 10 | parallel -j 8 "echo 4 4 | g_rms -f $DATA/ab_${iso}_${ratio}_{}_nosol_whole.xtc_c_fit.xtc -s ${iso}_${ratio}_nosol.tpr -o $output_dir/ab_${iso}_${ratio}_{}_nosol_whole_rmsd_backbone.xvg -noxvgr $TEST"
	
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
	seq 1 10 | parallel -j 8 "echo $calpha | g_rmsf -f $DATA/ab_${iso}_${ratio}_{}_nosol_whole.xtc_c_fit.xtc -s ${iso}_${ratio}_nosol.tpr -o $output_dir/${iso}_${ratio}_{}_rmsf.xvg -fit -res -ox $output_dir/ab_${iso}_${ratio}_{}_nosol_whole.xtc_c_fit -noxvgr $TEST" 
	clean "${iso}_${ratio}_rmsf"
}

mode='production'
TEST="-b 0"
if [ "$mode" == "production" ]; then
	echo "running in production mode"
	cd $PBS_O_WORKDIR
else
	echo "testing ..."
	#set externally bound variables
	ISO=scyllo
	RATIO=15
	ANALYSIS=dssp
	TEST="-b 1000 -e 1010"
fi

base_dir=`pwd`
DATA="/scratch/grace/inositol/abeta42/2/xtc"
SHM="/dev/shm/grace/"

echo `pwd`
${ANALYSIS} $ISO $RATIO $SHM

#remove everything in memory in case signals aren't trapped properly
rm -rf /dev/shm/*
