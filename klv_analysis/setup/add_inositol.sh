#!/bin/bash

# script used to generate a set of inositol + klvffae + disordered at 45:4 ratio
# Feb 16 2011

set -u
set -e

peptide="klv"
# peptide[1]="q8";
# peptide[2]="f8";

# size[0]=131;
# size[1]=145;
# size[2]=169;

for isomer in chiro; do
    for (( i=0 ; i < 25 ; i++ )); do 
		# mkdir $peptide/${isomer}/sys${i}
		# cp templates/${peptide}_ins_template.top $peptide/${isomer}/sys${i}/${peptide}_olig_${isomer}_${i}.top
        mkdir -p $isomer/sys${i}
		cp template.top $isomer/sys${i}/sys${i}.top
		genbox -cp ${peptide}_box/${peptide}_olig_${i}.gro -cs ~/gro/${isomer}_giant_box.gro -p $isomer/sys${i}/sys${i}.top -o $isomer/sys${i}/sys${i}.gro -maxsol 45
		genbox -cp $isomer/sys${i}/sys${i}.gro -cs ~/gro/tip3.gro -p $isomer/sys${i}/sys${i}.top -o $isomer/sys${i}/sys${i}_start.gro
		
	done
done
