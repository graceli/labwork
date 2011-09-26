#!/bin/sh
set -u
set -x
set -e

function preprocess_chain_hbonds {
	for iso in chiro scyllo water; do 
		for (( i=1; i<=10; i++ )); do
			for (( j=0; j<4; j++ )); do
				if [ -e "$iso/$i/chain_${j}_$((j+1))_hbonds.xvg" ]; then
					sed -e 's/@/#/g' $iso/$i/chain_${j}_$((j+1))_hbonds.xvg > ${iso}_sys${i}_chain_${j}_$((j+1))_hbonds.dat
				fi
			done
		done
	done
}

# unpack non chain_hbond files
# for analysis in rmsf; do 
#     for file in *15*rmsf*.tgz; do tar xvfz $file; done
# done

function unpack {
    for iso in scyllo chiro water; do
        for file in `ls *$iso*$ratio*$analysis.tgz`; do
            echo $file
            mkdir -p $analysis/$iso
            cp -p $file $analysis/$iso
            cd  $analysis/$iso
            tar xvfz $file
            cd -
        done
    done
}

# ratio=$1
# analysis=$2
# unpack

preprocess_chain_hbonds