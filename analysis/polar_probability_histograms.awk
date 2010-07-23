#/usr/bin/awk -f

# this script works for polar output files with 4 inositols
# klvffae specific
BEGIN{
	#initialize the histograms so don't have to do undefined values
	for(i=0;i<=10;i++){
		backbone_hist[i]=0;
		sidechain_lys_hist[i]=0;
		sidechain_glu_hist[i]=0;
	}	
}
{
	backbone_hist[$1]++;
	backbone_hist[$2]++;
	backbone_hist[$3]++;
	backbone_hist[$4]++;

	#lysines
	sidechain_lys_hist[$5]++;
	sidechain_lys_hist[$6]++;
	sidechain_lys_hist[$7]++;
	sidechain_lys_hist[$8]++;

	#glutamates
	sidechain_glu_hist[$9]++;
	sidechain_glu_hist[$10]++;
	sidechain_glu_hist[$11]++;
	sidechain_glu_hist[$12]++;
}
END{
	for(i=0;i<=5;i++){
		print i, backbone_hist[i]/NR/4, sidechain_lys_hist[i]/NR/4, sidechain_glu_hist[i]/NR/4,backbone_hist[i], sidechain_lys_hist[i], sidechain_glu_hist[i], NR;	
	}	
	
}

