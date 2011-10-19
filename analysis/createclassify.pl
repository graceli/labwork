#!/usr/bin/perl

print "usage: convert2gro.pl\n";

#@threads = (7,18,93,165,185,203);
my $num_threads=4;
my $num_segs=4;
my $residue="gly";
my $isomer="scy";

for(my $tnum=0; $tnum<=$num_threads; $tnum++){
	for(my $seg=0; $seg<=$num_segs; $seg++){
		#print "tpr/${residue}_${isomer}_${tnum}_${seg}.tpr\n";
		if(-e "tpr/${residue}_${isomer}_${tnum}_${seg}.tpr"){
			my $command = "dipep_class gro/${residue}_${isomer}_${tnum}_${seg}_part.gro rama/${residue}_${isomer}_${tnum}_${seg}_rama.xvg 2> /dev/null > analysis/${residue}_${isomer}_${tnum}_${seg}.class&";
			print $command,"\n";
			system($command);
		}
	}
}

