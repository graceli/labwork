#!/bin/sh

# analysis for Abeta17-42 protofibril for a single trajectory

proteingrp=1
insgrp=12
chain1=18
chain2=19
chain3=20
chain4=21
chain5=22
xvgr="-noxvgr"
start=1
end=20
#extract the number of nonpolar contacts between inositol and residues

function make_indices {
	echo -e "'SideChain'&aC*&!rACE\nsplitch17\nq" | make_ndx -f common/em.tpr -o common/ab_nonpolar.ndx 
	RETURN_CODE=$?
}

function nonpolar {
	for s in `seq $start $end`; do
		DATA=../sys${s}
		#mkdir sys${s}_o
		echo -e "'SideChain'&aC*&!rACE\nsplitch17\nq" | make_ndx -f $DATA/em.tpr -o $DATA/ab_nonpolar.ndx 
		for xtc in $DATA/sys${s}.xtc $DATA/sys${s}.part0002.xtc $DATA/sys${s}.part0003.xtc; do 
			echo "analyzing $xtc ..."
			for i in $chain1 $chain2 $chain3 $chain4 $chain5; do
				echo $i $insgrp | g_inositol_residue_nonpolar_v2 -f $xtc -s $DATA/em.tpr -n $DATA/ab_nonpolar.ndx -per_residue_contacts chain${i}_residue_np_contact.dat -per_inositol_contacts chain${i}_inositol_np_contact.dat
				 mv table.dat chain${i}_table.dat
			
				trap "exit 1" SIGTERM TERM KILL SIGINT
			done

			### evoke python script ###
			python abeta_analysis.py sys${s}.h5
			###########################
			rm *.dat
		done
	done
}

#TODO: hbonds

res_start=0
res_end=129
INS_grp=130
num=0
function hbonds {
	#for s in `seq $start $end`; do
	for s in `seq 1 20`; do
		DATA=../sys${s}
		for xtc in $DATA/sys${s}.xtc $DATA/sys${s}.part0002.xtc $DATA/sys${s}.part0003.xtc  $DATA/sys${s}.part0004.xtc; do
			for res_grp in `seq $res_start $res_end`; do
				echo $res_grp $INS_grp | g_hbond -f $xtc -s $DATA/em.tpr -n g_hbond_15.ndx -nonitacc -nomerge -num $res_grp $xvgr > /dev/null 2>&1 &
				let num=num+1
				if [ "$num" == "8" ]; then
					echo "waiting ..."
					wait
					num=0
				fi
				trap "exit 1" SIGTERM TERM KILL SIGINT
			done
			wait
#			python abeta_analysis.py sys${s}.h5
			python /home/grace/AnalysisScripts/abeta_analysis/abeta_analysis.py sys${s}.h5
		done
	done	
}



#make_indices
#nonpolar
cd $PBS_O_WORKDIR
hbonds
