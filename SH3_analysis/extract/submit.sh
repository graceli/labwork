#!/bin/sh

DATA=/scratch/p/pomes/grace/drSH3
num=0
for dir in `ls -d $DATA/PRIOR*`; do
	echo "trap 'rm -rf /dev/shm/*' INT TERM SIGINT; python $HOME/labwork/SH3_analysis/extract/extract_xtcs.py $dir" > extract_${num}.sh
	chmod u+x extract_${num}.sh
	let num=num+1
	# qsub -l nodes=1:compute-eth:ppn=8,walltime=10:00:00  extract_${num}.sh
done

for i in `seq 0 26`; do 
	qsub -l nodes=1:compute-eth:ppn=8,walltime=05:00:00  extract_${i}.sh
done

