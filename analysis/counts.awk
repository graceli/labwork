#!/usr/bin/awk -f
BEGIN{
	if(!COLUMN){
		print "undefined parameter variable COLUMN" > "/dev/stderr";
		exit;
	}	
	counts=0;
}
{
	if($COLUMN){
		counts++;
	}
}
END{
	print counts/NR, counts, NR;
}
