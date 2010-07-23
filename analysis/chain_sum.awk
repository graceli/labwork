#!/bin/awk -f
BEGIN{
	ubound=0;
	mono=0;
	two=0;
	three=0;
	four=0;
}
{
	ubound+=$1;
	one+=$2;
	for(i=4; i<=10;i++)
		mono_arr[i]+=$i;

	two+=$12;
	three+=$13;
	four+=$14;
}
END {
	user=ENVIRON["USER"];
	print "file analyzed :",user"@:"ENVIRON["PWD"]"/"FILENAME;
	print "Total number of inositols bound to n chains";
	print "n=0    n=1  n=2  n=3    n=4";
	print ubound, one, two, three, four;
	print "";
	print "breakdown of single chain:";
        for(i=4; i<=10;i++){
		print "#bound to",i-4,"groups=", mono_arr[i];
	}
			
}
