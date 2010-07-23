#!/usr/bin/awk -f
{
	#assumes that first column is time
	#collect data
	for(i=FIRST; i<=LAST; i++){
		line_total+=$i
	};
	print NR-1, line_total;
	line_total=0;
}
