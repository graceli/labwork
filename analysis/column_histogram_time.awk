#!/usr/bin/awk -f
BEGIN{
	#initialize histogram
	for(i=0;i<=20;i++){
		histogram[i]=0;
	}
	total_cols=0;
}
{
	#collect data
	for(i=2;i<=NF;i++){
		histogram[$i]++;
	};
	total_cols=NF-1;
}
END{
	#output histogram
	for(i=0;i<=20;i++){
		print i,histogram[i]/NR/total_cols, histogram[i], NR, total_cols
	}
}

