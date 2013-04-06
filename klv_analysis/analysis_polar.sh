#!/bin/sh
#PBS -l nodes=1:compute-eth:ppn=8,walltime=08:00:00,os=centos53computeA
#PBS -N analysis

set -u
set -e
set -x 
trap "clean; exit" INT TERM KILL

TEST='-dt 2'
CUTOFF=0.45
OUTPUT=/dev/null
NTASK=8

GROUP[15,0]=15
GROUP[15,1]=33
GROUP[45,0]=15 
GROUP[45,1]=63

#interpeptide nonpolar interactions
function pp_nonpolar {
	GRP=14
	task=0
	xtc=$1
	tpr=$2
	ndx=$3
	output_dir=$4/pp_nonpolar
	mkdir -p $output_dir
	for file in `ls $xtc/*.xtc`; do
		base=`basename $file .xtc`
        echo $GRP | g_pp_nonpolar -f $file -s $tpr -n $ndx -deffnm $output_dir/${base}_ -cutoff $CUTOFF $TEST > $OUTPUT 2>&1 &
        echo "starting process $task"
        let task=$task+1
        if [ "$task" == "$NTASK" ]; then
                wait
                task=0
        fi
	done
}

#inositol residue binding
function nonpolar_residue {
	GRP="14 12"
	task=0
	xtc=$1
	tpr=$2
	ndx=$3
	output_dir=$4/nonpolar_residue
	mkdir -p $output_dir
	for file in `ls $xtc/*.xtc`; do
		base=`basename $file .xtc`
		echo $GRP | g_inositol_residue_nonpolar_v2 -f $file -s $tpr \
			-n $ndx -per_residue_contacts $output_dir/${base}_per_residue_contact.dat \
			-per_inositol_contacts $output_dir/${file}_per_inositol_contacts.dat -dist $CUTOFF $TEST > $OUTPUT 2>&1 &
			
		let task=$task+1
		if [ "$task" == "$NTASK" ]; then
			wait
			task=0
		fi
	done
}

#polar analysis
#peptide - peptide 
#peptide - inositol
function polar_residue {
	task=0
	xtc=$1
	tpr=$2
	ndx=$3
	output_dir=$4/polar
	NPEP=$5
	NINOS=$6
	GRP1=$7
	GRP2=$8
	mkdir -p $output_dir
	for file in `ls $xtc/*.xtc`; do
		base=`basename $file .xtc`
		seq $GRP1 $GRP2 | g_parse_index_oct21 -f $xtc/$base -s $tpr -n $ndx -num_peptides $NPEP -num_inositol $NINOS -deffnm $output_dir/${base}_ $TEST 2> $OUTPUT >&2 &
			
		let task=$task+1
		if [ "$task" == "$NTASK" ]; then
			wait
			task=0
		fi
	done
}

#cluster size calculation with g_clustsize
function cluster {
	GRP=1 #protein
	task=0
	xtc=$1
	tpr=$2
	ndx=$3
	output_dir=$4/cluster
	mkdir -p $output_dir
	for file in `ls $xtc/*.xtc`; do
		base=`basename $file .xtc`
	    echo $GRP | g_clustsize -f $file -s $tpr -n $ndx -nc $output_dir/${base}_nclust.xvg -cut $CUTOFF -noxvgr $TEST > $OUTPUT 2>&1 &
		let task=$task+1
		if [ "$task" == "$NTASK" ]; then
			wait
			task=0
		fi
	done
}

#secondary structure calculation with DSSP
function dssp {
	GRP=1 #protein
	xtc=$1
	tpr=$2
	ndx=$3
	output_dir=$4/dssp
	mkdir -p $output_dir
	for file in `ls $xtc/*.xtc`; do
		base=`basename $file .xtc`
		echo $GRP | do_dssp -f $file -s $tpr -o $output_dir/${base}_o -sc $output_dir/${base}_sc $TEST > $OUTPUT 2>&1 &
	done	
}

function clean {
#	rm -rf dssp cluster polar nonpolar_residue pp_nonpolar
	cd /dev/shm
	tar cvfz "analysis${1}.tgz" analysis
	cp analysis*.tgz $base_dir
	rm -rf analysis analysis.tgz
	cd -
}

#do all analysis for systems in the presence of inositol
function do_inositol {	
	system=$1
	for sys in $system; do
		xtc="${sys}to4/xtc/0-200ns"
		tpr="common/${sys}to4_nosol.tpr"
		ndx="common/${sys}to4_analysis.ndx"
		output_base="/dev/shm/analysis"

		NPEP=4; NINOS=$sys
		polar_residue $xtc $tpr $ndx $output_base $NPEP $NINOS ${GROUP[$sys,0]} ${GROUP[$sys,1]}
	done
}

#do all analysis for systems in the absence of inositol
function do_with_water {
	xtc=...
	TPR=

	pp_nonpolar $xtc
	wait
	nonpolar_residue $xtc
	wait
	cluster $xtc
	wait
	dssp $xtc
	wait	
}

# cd $PBS_O_WORKDIR
SYSTEM=45
base_dir=`pwd`
do_inositol $SYSTEM
rm /dev/shm/analysis/*/\#*
rm \#*
clean ${SYSTEM}_polar

