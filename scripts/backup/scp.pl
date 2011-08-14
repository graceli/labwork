#!/usr/bin/perl
if(@ARGV <1){
	print "perl scp.pl <dir>";
	exit;
}

$curdir=$ARGV[0];
#@threads = (7,18,93,165,185,203);
#foreach(@threads){
for(my $num=0; $num<=3; $num++){
	my $scp_xtc = "scp ccb:/projects/pomes/grace/inositol/gly_dipep/epi/inos_${_}${num}/${num}.xtc $curdir/ala_scy_${_}${num}.xtc";
	my $scp_tpr = "scp ccb:/projects/pomes/grace/inositol/gly_dipep/epi/inos_${_}${num}/${num}.tpr $curdir/ala_scy_${_}${num}.tpr";
	#print $scp_xtc,"\n",$scp_tpr,"\n";
	system($scp_xtc);
	system($scp_tpr);
}
#}
