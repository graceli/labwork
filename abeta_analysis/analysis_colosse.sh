#!/bin/sh
#PBS -l nodes=1:compute-eth:ppn=8,walltime=02:00:00,os=centos53computeA
#PBS -N analysis

#set -u
#set -e
set -x

trap 'exit' TERM INT SIGINT EXIT SIGKILL SIGSTOP SIGTERM

# res_start=0
# res_end=129
# INS_grp=130
# num=0
# function hbonds {
# 	iso=$1
#     ratio=$2
# 	output_dir=$3/hbonds
# 	mkdir -p $output_dir
# 		
# 	for s in `seq 1 10`; do
# 		xtc="ab_${iso}_${ratio}_${s}_nosol_whole.xtc_c_fit"
# 		if [ -e "$DATA/${xtc}.xtc" ]; then
# 			mkdir -p $output_dir/$s
# 			seq $res_start $res_end | parallel -j 8 "echo {} $INS_grp | g_hbond -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n g_hbond_${ratio}_${iso}.ndx -nonitacc -nomerge -num $output_dir/$s/{} $xvgr $TEST > /dev/null 2>&1"
# 		fi
# 		# python /home/grace/AnalysisScripts/abeta_analysis/abeta_analysis.py sys${s}.h5
# 	done
# 	clean "${iso}_${ratio}_hbonds"
# }

# proteingrp=1
# insgrp=12
# chain1=17
# chain2=18
# chain3=19
# chain4=20
# chain5=21
# xvgr="-noxvgr"
# # start=1
# # end=20
# function nonpolar {
# 	iso=$1
#     ratio=$2
# 	output_dir=$3/nonpolar
# 	mkdir -p $output_dir
# 	
# 	 echo -e "'SideChain'&aC*&!rACE\nsplitch16\nq" | make_ndx -f ${iso}_${ratio}_nosol.tpr -o ab_${ratio}_nonpolar.ndx
# 	
# 	for s in `seq 1 10`; do
# 		xtc="ab_${iso}_${ratio}_${s}_nosol_whole.xtc_c_fit"
# 		if [ -e "$DATA/${xtc}.xtc" ]; then
# 			seq $chain1 $chain5 | parallel -j 5 "echo {} $insgrp | g_inositol_residue_nonpolar_v2 -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n ab_${ratio}_nonpolar.ndx -per_residue_contacts $output_dir/chain{}_residue_np_contact.dat -per_inositol_contacts $output_dir/chain{}_inositol_np_contact.dat -per_residue_table chain{}_table.dat $TEST"
# 		fi
# 	done
# 	clean "${iso}_${ratio}_nonpolar"
# }


function dssp {
	echo 1 | do_dssp -f $NAME -s $TPR -o ${NAME}_ss -sc ${NAME}_sc
}

chain_start=0
chain_end=3
function chain_hbonds {
	xtc="ab_${iso}_${ratio}_${s}_nosol_whole.xtc_c_fit"
	for ch in `{$chain_start..$chain_end}`; do
		let next=ch+1
		echo $ch $next | g_hbond -f $NAME -s $TPR -n $DATA/chain.ndx -nonitacc -nomerge -num chain_${ch}_${next}_hbonds $TEST
	done
}

# calculate the rmsf
function rmsf {
    # backbone fitting for specific parts of the peptide
	echo "C-alpha" | g_rmsf -f $NAME -s $TPR -o sys${SYS}_${TAG}_rmsf.xvg -fit -res -ox sys${SYS}_${TAG}_ox.xvg -noxvgr $TEST 
}

# calculate the rmsd of the protein using the nmr structure as a reference
function rmsd {
	echo 1 1 | g_rms  -o sys${SYS}_whole_rmsd_protein.xvg -noxvgr $TEST
}

base_dir=`pwd`
DATA="/rap/uix-840-ac/grace/abeta/42/glucose/xtc"

# SYS=$SYSN
# ANALYSIS=$ANALY
SYS=1
ANALYSIS=dssp
TEST="-dt 1000"
TAG="nosol_whole"
NAME="$DATA/sys${SYS}_${tag}.xtc"
TPR="$DATA/protein_sugar.tpr"
echo "For system $SYS"
echo "Performing analysis $ANALYSIS"

mkdir analysis/$SYS; cd analysis/$SYS
mkdir $ANALYSIS; cd $ANALYSIS
${ANALYSIS} $SYS

