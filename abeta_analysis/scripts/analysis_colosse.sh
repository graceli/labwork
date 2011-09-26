#!/bin/bash
#$ -N ab_analysis
#$ -P uix-840-ac
#$ -A uix-840-ac
#$ -l h_rt=5:00:00
#$ -pe default 16
#$ -q med
#$ -S /bin/bash
#$ -cwd
#$ -notify

#module load tools/gnu-parallel/20110522

set -u
set -e
set -x

trap 'exit 1' TERM INT SIGINT EXIT SIGKILL SIGSTOP SIGTERM

function dssp {
	echo 1 | do_dssp -f $NAME -s $TPR -o ${NAME}_ss -sc ${NAME}_sc
}

chain_start=0
chain_end=4
function chain_hbonds {
	for (( i=${chain_start}; i < ${chain_end}; i++ )); do
		let next=i+1
		echo $i $next | g_hbond -f $DATA/$NAME -s $TPR -n $DATA/chain.ndx -nonitacc -nomerge -num ${NAME}_chain_${i}_${next}_hbonds $TEST
	done
}

# calculate the rmsf
function rmsf {
    	# backbone fitting for specific parts of the peptide
	echo "C-alpha" | g_rmsf -f $DATA/$NAME -s $TPR -o ${NAME}_rmsf.xvg -fit -res -ox ${NAME}_ox.xvg -noxvgr $TEST 
}

# calculate the rmsd of the protein using the nmr structure as a reference
function rmsd {
	echo 1 1 | g_rms -f $DATA/$NAME -s $TPR -o ${NAME}_whole_rmsd_protein.xvg -noxvgr $TEST
}

function run_analysis {
	analysis=$1	
	TEST="-b 0"
	for ((s=1; s<=10; s++)); do
		TAG="nosol_whole"
		NAME="sys${s}_${TAG}"
		TPR="$DATA/protein_sugar.tpr"
		${analysis}
	done
}

base_dir=`pwd`
DATA="/rap/uix-840-ac/grace/abeta/42/glucose/xtc"

echo "in $PWD"

for ANALYSIS in rmsd; do
	echo "Performing analysis $ANALYSIS"

	if [ ! -e "$ANALYSIS" ]; then
		echo "$ANALYSIS does not exist ... making directory"
		mkdir $ANALYSIS
	fi

	cd $ANALYSIS
	run_analysis $ANALYSIS &
	cd ..
done

