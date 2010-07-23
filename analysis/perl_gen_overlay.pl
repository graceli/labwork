#!/usr/bin/perl
#

use strict;

my @backbone=qw(CO0CO2 NH1CO2 CO0NH3 NH1NH3); #backbone combos
my $bb; 
my $dipep = "ala";  #dipeptide short name
my $isomer = "scy"; #made up inositol isomer short name
my $short = substr($dipep,0,1); #first character of dipeptide name
foreach $bb (@backbone){
	for(my $i=12; $i<=14; $i++){
		#generate gnuplot 4.2 plot statements
		if(-e "${dipep}_${isomer}_${bb}_${i}.dat" && -s "${dipep}_${isomer}_${bb}_${i}.dat") {
			print "plot [-180:180][-180:180] \'${dipep}_${isomer}_${bb}_${i}.dat\' ";
			#if(-s "${dipep}_${isomer}_${bb}_${i}.dat" > 5000) {
			#	print "w dots"; 
			#}else{
			#	print "w points";
			#}
 			print ", \'${short}wat_contour.dat\' w l \n";
		} else {
			print "unset grid\n";
			print "plot [0:0.1][0:0.1] 0 w dots\n";
			print "set grid\n";
		}
	} 
}
