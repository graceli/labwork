#!/usr/bin/awk -f
{
	#collect data
	for(i=1;i<=NF;i++){
		line_total+=$i
	};
	print NR-1, line_total;
	line_total=0;
}

