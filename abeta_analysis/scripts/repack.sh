#!/bin/sh

set -u 
set -x

function repack {
for ratio in 15 64; do
	for iso in scyllo chiro glycerol glucose water; do
		for analysis in rmsf rmsd nonpolar hbonds chain_hbonds; do
			FILE="analysis_${iso}_${ratio}_${analysis}.tgz"
			DIR=$ratio/$iso/$analysis
			echo "Extracting $FILE ... into $DIR"
			mkdir -p $DIR
			cp $FILE $DIR
			cd $DIR
			tar xvfz $FILE 
			mv analysis/$analysis/* . 
			rm -rf analysis 
			cd ../../../
		done
	done
done
}


function repack2 {
for ratio in 15 64; do
	for iso in scyllo chiro glycerol glucose water; do
		for analysis in chain_hbonds; do
			FILE="analysis_${iso}_${ratio}_${analysis}.tgz"
			DIR=$ratio/$iso/$analysis
			echo "Extracting $FILE ... into $DIR"
			mkdir -p $DIR
			cp $FILE $DIR
			cd $DIR
			tar xvfz $FILE 
			mv analysis/hbonds/* . 
			rm -rf analysis 
			rm -f $FILE
			cd ../../../
		done
	done
done
}

repack2
