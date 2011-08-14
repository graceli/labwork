#!/usr/bin/perl

$curdir=$ARGV[0];
@threads = (7,18,93,165,185,203);
foreach(@threads){
	for(my $num=0; $num<=3; $num++){
		my $command = "g_rama -f ala_scy_${_}${num}.xtc -s ala_scy_${_}${num}.tpr -o ala_scy_${_}${num}_rama.xvg";
		#print $command,"\n";
		system($command);
	}
}
