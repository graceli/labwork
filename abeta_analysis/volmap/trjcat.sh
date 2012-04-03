#!/bin/sh
#PBS -l nodes=1:ppn=8,walltime=10:00:00,os=centos53computeA
#PBS -N volmap

cd $PBS_O_WORKDIR

DATA=../xtc
FIT=0
OUTPUT=1

for ratio in 15 64; do 
	for iso in glycerol chiro scyllo; do
		xtc_file="${iso}_${ratio}_volmap_all"
		#for i in `seq 1 100`; do echo "c"; done | trjcat -f $DATA/*${iso}*${ratio}*c_fit.xtc -o $xtc_file -cat -keeplast -settime
		echo $FIT $OUTPUT | trjconv -f $xtc_file -s ../${iso}_${ratio}_nosol.tpr -fit rot+trans -o ${xtc_file}_fit -n ../${iso}_${ratio}_nosol.ndx
		if [ "$?" == "0" ]; then
			#rm $xtc_file
			echo "finished fitting on $xtc_file"
		else
			echo "Error trjconving $xtc_file"
		fi
	done
done
