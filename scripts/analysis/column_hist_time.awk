#!/usr/bin/awk -f 

#NCOLS -- number of columns in the file -- is an input parameter

BEGIN {
	for(b=0;b<=20;b++){
		hist[b]=0;
	}
}
{
	#first column is time
	for(i=2;i<=NCOLS+1;i++){
		hist[$i]++;	
	}
}
END {
	for(b=0;b<=20;b++){
		print b, hist[b]/NR/NCOLS,hist[b];
	}
}

