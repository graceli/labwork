#!/usr/bin/awk -f
BEGIN{
	#initialize histogram
	for(i=0;i<=20;i++){
		histogram[i]=0;
	}
}
{
	#collect data
	for(i=2;i<=NF;i++){
		line_total+=$i
	};
	histogram[line_total]++;
	line_total=0;
}
END{
	#output histogram
	for(i=0;i<=20;i++){
		print i,histogram[i]/NR, histogram[i], NR
	}
}

