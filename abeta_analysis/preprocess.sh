#!/bin/sh

function sort_data {
	for iso in glycerol; do
		for r in 15 64; do 
			if [ ! -e "${iso}_${r}" ]; then
				mkdir ${iso}_${r}
			fi
			mv data/*${iso}_${r}* ${iso}_${r}
			# cp config.ini ${iso}_${r}
		done
	done
}

function read_into_pytables {
	for iso in glycerol; do
		for r in 15 64; do
			# other analysis nonpolar hbond
	        for analysis in rmsf chain_hbonds; do
				cd ${iso}_${r}
				tar xvfz analysis_${iso}_${r}_${analysis}.tgz
				echo "moving stuff in analysis/${analysis} to $PWD"
	            mv analysis/${analysis} .
	        	echo "saving into pytables ..."
	        	# python ~/AnalysisScripts/abeta_analysis/abeta_analysis.py ../${iso}_${r}.h5
	        	# rm *.dat *.xvg
	            #rm -rf analysis
				cd ../
	        done
		done
	done
}

# sort_data
read_into_pytables