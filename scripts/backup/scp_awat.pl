#!/usr/bin/perl
@t = (185, 203, 7);
@l = (4, 4, 5);
for($k =0; $k<3; $k++){
	for($i=0; $i<$l[$k]; $i++){
		$scp="scp ccb:/projects/pomes/grace/inositol/ala_dipep/water/water_dipep_threads/awat$t[$k]_${i}/${i}.tpr awat$t[$k]_${i}.tpr";
		print $scp,"\n";
		if($ARGV[0] eq "-exec") {
			system($scp);
		}
	}
}
