#/usr/bin/awk -f

# this script works for polar output files with 4 inositols
# klvffae specific

#LEN -- number of residues per peptides
#NRES -- number of peptides
BEGIN{
	if(!LEN || !NRES){
		print "collapse.awk: parameters LEN and NRES are not set">"/dev/stderr";
		exit;
	}
}
{
	for(n=1; n<=LEN; n++){
		for(i=(n-1)*NRES+2;i<=n*NRES+1;i++){
			hist[n]+=$i;
		}
	}
}
END{
	for(i=1;i<=LEN;i++){
		print i, hist[i]/NR, hist[i], NR;
	}
}

