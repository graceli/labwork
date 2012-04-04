#!/bin/bash
#$ -N ab_analysis
#$ -P uix-840-ac
#$ -A uix-840-ac
#$ -l h_rt=10:00:00
#$ -pe default 16
#$ -q med
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
set -x

trap "clean; exit $?" TERM INT SIGINT EXIT SIGKILL SIGSTOP SIGTERM

function clean {
	cd /dev/shm/grace
	tar cvfz ${base_dir}/analysis.tgz *
	rm -rf /dev/shm/grace
}
#centers and puts entire fibril in box
function preprocess {
	for i in `seq 5 9`; do
		echo 2 1 | trjconv -f sys$i/sys${i}_prod.xtc -s sys$i/sys${i}_prod.tpr -n $DATA/protein_sugar.ndx -o $DATA/sys${i}_c_res_nosol -center -pbc res
	done
}

function dssp {
	# break dssp analysis to write to different directories because temp file names can conflict
	# todo: fix these bd temp files??
	# note: xvgr is useful for dssp analysis

	mkdir /dev/shm/grace
	cd /dev/shm/grace

	export DSSP=/home/grace/labwork/gromacs/dssp_ana/dsspcmbi

	if [ ! -e "$1" ]; then
		mkdir $1
	fi
	cd $1
	echo 1 | do_dssp -f $DATA/$NAME -s $DATA/$TPR -o ${NAME}_ss -sc ${NAME}_sc $TEST -dt 10 &
	cd ../
}

#computes number of interchain hydrogen bonds (between adjacent chains)
chain_start=0
chain_end=4
function chain_hbonds {
	for (( i=${chain_start}; i < ${chain_end}; i++ )); do
		let next=i+1
		echo $i $next | g_hbond -f $DATA/$NAME -s $DATA/$TPR -n $DATA/chain.ndx -nonitacc -nomerge -num ${NAME}_chain_${i}_${next}_hbonds -noxvgr $TEST &
	done
	wait
}

# calculate the rmsf for each residue in fibrl
function rmsf {
    	# backbone fitting for specific parts of the peptide
	echo "C-alpha" | g_rmsf -f $DATA/$NAME -s $DATA/$TPR -o ${NAME}_rmsf.xvg -fit -res -ox ${NAME}_ox.xvg -noxvgr $TEST &
}

# calculate the rmsd of the protein using the nmr structure as a reference
function rmsd {
	echo 1 1 | g_rms -f $DATA/$NAME -s $DATA/$TPR -o ${NAME}_whole_rmsd_protein.xvg -fit rot+trans -noxvgr $TEST &
}

function gyration {
	echo 1 | g_gyrate -f $DATA/$NAME -s $DATA/$TPR -o ${NAME}_rg_protein.xvg -noxvgr $TEST 2> rg.out >&2 &
} 

# task function to run analysis
function run_analysis {
	analysis=$1	
	for ((s=1; s<=10; s++)); do
		NAME="${s}_${TAG}"
		TPR="protein_glca.tpr"
		
		if [ -e "$DATA/${NAME}.xtc" ]; then
			echo "running task $s"
			${analysis} ${s}
		else	
			echo "trajectory ${NAME}.xtc does not exist"
	        fi
	done
}

# run all analysis at once
function batch_run {
	for ANALYSIS in dssp; do
		echo "Performing analysis $ANALYSIS"

		if [ ! -e "$ANALYSIS" ]; then
			echo "$ANALYSIS does not exist ... making directory"
			mkdir $ANALYSIS
		fi

		cd $ANALYSIS
		run_analysis $ANALYSIS
		wait
		cd ..
	done
}

base_dir=`pwd`
DATA="/rap/uix-840-ac/grace/abeta/42/glucose_Protein_GLCA"
TAG="final"
TEST="-b 0"

echo "in $PWD"
echo "running app"

#exec app from command line
batch_run

