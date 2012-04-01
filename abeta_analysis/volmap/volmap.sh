#!/bin/sh

DATA=../xtc
for ratio in 15 64; do 
	for iso in glycerol chiro scyllo; do
		xtc_file="${iso}_${ratio}_volmap_all.xtc"
		for i in `seq 1 100`; do echo "c"; done | trjcat -f $DATA/*${iso}*${ratio}*c_fit.xtc -o $xtc_file -cat -keeplast -settime
		
		is_cer=0
		if [ "$iso" == "glycerol" ]; then
			is_cer="1"
		fi

		qsub -v XTC="$xtc_file",GRO="${iso}_${ratio}_nosol.gro",is_cer=$is_cer volmap_template.sh
	done
done
