#!/usr/bin/awk -f

## col1  col2  col3             col4              col5
## intra inter normalized intra normalized inter  total number of frames

{
	if($2 ~ /[0-9]+/){
		#intra-HB
		total_intra+=$2;

		#inter-HB
		total_inter+=$3;

		#inositol-peptide
		total_inos+=$4;
	}	
}

END{
#	print "total_inter = ", total_inter;
#	print "total_inter(per ps per chain) = ", total_inter/(NR*4);
#	print "total_intra = ", total_intra;
#	print "total_intra(per ps per chain) = ", total_intra/(NR*4);

	norm_inter= total_inter/(NR*c_npep);
	norm_intra = total_intra/(NR*c_npep);
	norm_inos = total_inos/(NR*nins);
	print total_intra,total_inter, total_inos, norm_intra, norm_inter, norm_inos, NR; # NR = number of lines in a file read
}

