#!/usr/bin/perl

@xtc = <xtc/*.xtc>;
@tpr = <tpr/*.tpr>;

for(my $i=0; $i<@xtc; $i++){
	$command="g_rama -f $xtc[$i] -s $tpr[$i] -o $xtc[$i] -noxvgr";
	print $command,"\n";
	if($ARGV[0] eq "-exec"){
		system($command);
	}

}
