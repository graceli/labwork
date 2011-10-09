#!/bin/bash
#$ -N ab_analysis
#$ -P uix-840-ac
#$ -A uix-840-ac
##$ -l h_rt=00:15:00
##$ -pe default 16
##$ -q med
#$ -S /bin/bash
#$ -cwd
#$ -notify

#module load tools/gnu-parallel/20110522
module load compilers/intel/11.1.059
module load mpi/openmpi/1.3.4_intel
export OMP_NUM_THREADS=$NSLOTS
# echo "Got $NSLOTS processors."
. /home/grace/.gmx

set -u
set -e
set -x

trap "exit $?" TERM INT SIGINT EXIT SIGKILL SIGSTOP SIGTERM


function dssp {
	# break dssp analysis to write to different directories because temp file names can conflict
	# todo: fix these bad temp files??
	export DSSP=/home/grace/labwork/gromacs/dssp_ana/dsspcmbi

	if [ ! -e "$1" ]; then
		mkdir $1
	fi
	cd $1
	echo 1 | do_dssp -f $DATA/$NAME -s $DATA/$TPR -o ${NAME}_ss -sc ${NAME}_sc &
	cd ../
}

chain_start=0
chain_end=4
function chain_hbonds {
	for (( i=${chain_start}; i < ${chain_end}; i++ )); do
		let next=i+1
		echo $i $next | g_hbond -f $DATA/$NAME -s $DATA/$TPR -n $DATA/chain.ndx -nonitacc -nomerge -num ${NAME}_chain_${i}_${next}_hbonds $TEST
	done
}

# calculate the rmsf
function rmsf {
    	# backbone fitting for specific parts of the peptide
	echo "C-alpha" | g_rmsf -f $DATA/$NAME -s $DATA/$TPR -o ${NAME}_rmsf.xvg -fit -res -ox ${NAME}_ox.xvg -noxvgr $TEST 
}

# calculate the rmsd of the protein using the nmr structure as a reference
function rmsd {
	echo 1 1 | g_rms -f $DATA/$NAME -s $DATA/$TPR -o ${NAME}_whole_rmsd_protein.xvg -noxvgr $TEST &
}

function run_analysis {
	analysis=$1	
	for ((s=1; s<=10; s++)); do
		NAME="sys${s}_${TAG}"
		TPR="protein_sugar.tpr"
		${analysis} ${s}
		echo "running task $s"
	done
}

base_dir=`pwd`
DATA="/rap/uix-840-ac/grace/abeta/42/glucose/xtc"
TAG="nosol_whole_res"
TEST=""

echo "in $PWD"

for ANALYSIS in rmsf rmsd; do
	echo "Performing analysis $ANALYSIS"

	if [ ! -e "$ANALYSIS" ]; then
		echo "$ANALYSIS does not exist ... making directory"
		mkdir $ANALYSIS
	fi

	cd $ANALYSIS
	run_analysis $ANALYSIS
	echo "waiting"
	wait
	cd ..
done

