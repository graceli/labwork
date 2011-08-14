#!/usr/bin/awk -f
# this script computes the probability distribution of inter-HB, and intra-HB (averaged over time)
# and computes separate distributions for the first half and the second half of a trajectory
# c_npep is set by -v c_npep=X on the command line

{
	#build the histograms
	if(!/^[a-zA-Z]+/) {
		intra_hist[$2]++;
		inter_hist[$3]++;
	}
}

END {
	# output the histograms normalized by the total number of snapshots	
	for(i=0;i<=24;i++){
		printf "%d %f\n",i,inter_hist[i]/(NR*c_npep) > "inter_dist_pd.hist";
	}

	for(i=0;i<=12;i++){
		printf "%d %f\n",i,intra_hist[i]/(NR*c_npep) > "intra_dist_pd.hist";
	}
	##print NR;
}
