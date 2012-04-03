#!/bin/sh
#PBS -l nodes=1:ppn=8,walltime=10:00:00
#PBS -N analysis

. ~/.gmx_407

set -u
set -x
#set -e

#extract the number of nonpolar contacts between inositol and residues
# function make_indices {
# 	echo -e "'SideChain'&aC*&!rACE\nsplitch17\nq" | make_ndx -f common/em.tpr -o common/ab_nonpolar.ndx 
# 	RETURN_CODE=$?
# }

function clean {
    cd /dev/shm/grace
	tar cvfz analysis_${1}.tgz *
	cp analysis_${1}.tgz $base_dir
	rm -rf /dev/shm/*
}


# Calculates the number of nonpolar contacts of solute to residue for each chain
# In the aggregate.  Don't really remember why I did it like this
proteingrp=1
# insgrp can be glycerol or glucose
insgrp=5
chain1=0
chain5=4
xvgr="-noxvgr"
# start=1
# end=20
function nonpolar {
	iso=$1
    ratio=$2
	output_dir=$3/nonpolar
	mkdir -p $output_dir
	
    # We don't need make the index over and over again for the peptides, because they are all the same
	 #echo -e "'SideChain'&aC*&!rACE\nsplitch16\nq" | make_ndx -f ${iso}_${ratio}_nosol.tpr -o ab_${ratio}_nonpolar.ndx
	
	for s in `seq 0 9`; do
		xtc="${s}_final"
		if [ -e "$DATA/${xtc}.xtc" ]; then
			seq $chain1 $chain5 | parallel -j 5 "echo {} $insgrp | g_inositol_residue_nonpolar_v2 -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n ab_${ratio}_${iso}_nonpolar.ndx -per_residue_contacts $output_dir/${s}_chain{}_residue_np_contact.dat -per_inositol_contacts $output_dir/${s}_chain{}_inositol_np_contact.dat -per_residue_table $output_dir/${s}_chain{}_table.dat -per_inositol_phe_contacts $output_dir/per_inositol_phe_contacts.dat -FF_info $output_dir/ff_vs_t.dat -com_dist_xvg $output_dir/per_inositol_phe_com_dists.dat $TEST"
		fi
	done
    #clean "${iso}_${ratio}_nonpolar"
}

# calculate the rmsd of the protein using the nmr structure as a reference
function rmsd {
    iso=$1
    ratio=$2
    output_dir=$3/rmsd
    mkdir -p $output_dir

    seq 0 9 | parallel -j 8 "echo 1 1 | g_rms -f $DATA/{}_final.xtc -s ${iso}_${ratio}_nosol.tpr -o $output_dir/{}_rmsd.xvg -noxvgr $TEST"

    clean "${iso}_${ratio}_rmsd"
}

# Calculates the number of hydrogen bonds made with each residue.
# To calculate this I run g_hbond each time for each residue in the pentamer
# This is a bit of a hack to get g_hbond to work.
res_start=0
res_end=129
INS_grp=130
num=0
function hbonds {
    iso=$1
    ratio=$2
    output_dir=$3/hbonds
    mkdir -p $output_dir
  
    for s in `seq 0 9`; do
      xtc="${s}_final"
      if [ -e "$DATA/${xtc}.xtc" ]; then
          mkdir -p $output_dir/$s
          seq $res_start $res_end | parallel -j 8 "echo {} $INS_grp | g_hbond -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n g_hbond_${ratio}_${iso}.ndx -nonitacc -nomerge -num $output_dir/$s/{} $xvgr $TEST > /dev/null 2>&1"
      fi
      # python /home/grace/AnalysisScripts/abeta_analysis/abeta_analysis.py sys${s}.h5
    done

    # Cleaning up
    clean "${iso}_${ratio}_hbonds"
}

chain_start=0
chain_end=3
function chain_hbonds {
    iso=$1
    ratio=$2
    output_dir=$3/chain_hbonds
    mkdir -p $output_dir
  
    for s in `seq 0 9`; do
      xtc="${s}_final"
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

function dssp {
    export DSSP=/home/p/pomes/grace/src/dssp_ana/dsspcmbi

    for i in `seq 0 9`; do 
      if [ ! -e "$SHM/$i" ]; then
          mkdir -p $SHM/dssp/$i
      fi
    done

    seq 0 9 | parallel -j 8 "cd $SHM/dssp/{}; echo 1 | do_dssp -f $DATA/{}_final -s $base_dir/${ISO}_${RATIO}_nosol.tpr -o ab_${ISO}_${RATIO}_{}_ss -sc ab_${ISO}_${RATIO}_{}_sc $TEST 2>&1"

    clean "${ISO}_${RATIO}_dssp"
}

# calculate the rmsf
function rmsf_calpha {
    iso=$1
    ratio=$2
    output_dir=$3/rmsf_calpha
    mkdir -p $output_dir
    calpha=3

    # backbone fitting for specific parts of the peptide
    seq 0 9 | parallel -j 8 "echo $calpha | g_rmsf -f $DATA/{}_final -s ${iso}_${ratio}_nosol.tpr -o $output_dir/${iso}_${ratio}_{}_rmsf.xvg -fit -res -ox $output_dir/{}_final -noxvgr $TEST" 
    
    clean "${iso}_${ratio}_rmsf"
}

source ~/.gmx_407

mode='test'

TEST="-b 0"
if [ "$mode" == "production" ]; then
	echo "running in production mode"
	cd $PBS_O_WORKDIR
else
	echo "testing ..."
	#set externally bound variables
	ISO=glycerol
	RATIO=15
	ANALYSIS=nonpolar
	TEST="-b 1000 -e 1010"
fi

#trap 'rm -rf /dev/shm/*; echo "last process ended with retcode=$?"' TERM INT SIGINT EXIT SIGKILL SIGSTOP SIGTERM

base_dir=`pwd`
DATA="$SCRATCH/inositol/abeta42/current/$RATIO/${ISO}_nonsolvent"
SHM="/dev/shm/grace/"

mkdir $SHM

# run the analysis
${ANALYSIS} $ISO $RATIO $SHM
