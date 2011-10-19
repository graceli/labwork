#!/usr/bin/perl
@t = (93, 165, 185, 203, 7);
@l = (4,4,10,10,10);
for($k =0; $k<5; $k++){
	for($i=0; $i<$l[$k]; $i++){
		$tpr="rsync --rsh=ssh --rsync-path=/tools/rsync/2.6.9/bin/rsync -av --progress ccb:/projects/pomes/grace/inositol/ala_dipep/water/water_dipep_threads/awat$t[$k]_${i}/${i}.tpr tpr/awat$t[$k]_${i}.tpr";
		$xtc="rsync --rsh=ssh --rsync-path=/tools/rsync/2.6.9/bin/rsync -av --remove-sent-files --progress ccb:/projects/pomes/grace/inositol/ala_dipep/water/water_dipep_threads/awat$t[$k]_${i}/${i}.xtc xtc/awat$t[$k]_${i}.xtc";
		print STDERR "$tpr\n$xtc\n";

		if($ARGV[0] eq "-exec") {
			system($tpr);
			system($xtc);
		}
	}
}
