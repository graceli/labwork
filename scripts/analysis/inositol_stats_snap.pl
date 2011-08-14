#!/usr/bin/perl
###  This version of inositol_stats prints out the ramachandran maps for all the bidentate 
#states with the file and frame number of the snapshot for ease of locating the snapshot in a trajectory
#
#
if($ARGV[0] eq "-h"){
	print "usage: perl inositol_stats.pl \n";
	print "Collect and calculate equilibrium data in form of *.class files in the current directory\n";
	print "Characterization of binding equilibrium: % bidentate, mono, and none-binding populations, peptide inositol network, and multibound populations statistics are outputted\n"; 
	print "Dihedral (phi,psi) are sorted into different data files according to their binding state\n";
	exit;
}

@ARGV = <*.class>;

my $base=substr($ARGV[0],0,7);

open(BID1,">${base}_CO0CO2_12_snap.dat");
open(BID2,">${base}_CO0NH3_12_snap.dat");
open(BID3,">${base}_NH1CO2_12_snap.dat");
open(BID4,">${base}_NH1NH3_12_snap.dat");

open(BID5,">${base}_CO0CO2_13_snap.dat");
open(BID6,">${base}_CO0NH3_13_snap.dat");
open(BID7,">${base}_NH1CO2_13_snap.dat");
open(BID8,">${base}_NH1NH3_13_snap.dat");

open(BID9,">${base}_CO0CO2_14_snap.dat");
open(BID10,">${base}_CO0NH3_14_snap.dat");
open(BID11,">${base}_NH1CO2_14_snap.dat");
open(BID12,">${base}_NH1NH3_14_snap.dat");

print "gathering done on the following files:\n @ARGV\n\n";
my @avg_all_files=(0,0,0,0,0,0,0);
my @stats_names=qw(bidentate  mono  none  diid  multi_bound bound_stoic kd);
foreach $filename (@ARGV){
	print "### $filename ###\n";
	my @stats = gather($filename);
	for(my $i=0; $i<@stats; $i++){
		$avg_all_files[$i]+=$stats[$i];
		printf '%s=%.3f%s',$stats_names[$i],$stats[$i],"\n";
	}
	print "\n";
}
print "### summary ###\n";
for(my $i=0; $i<@avg_all_files-1; $i++){
	my $val = $avg_all_files[$i];
	my $avg = $val/@ARGV;
	printf '%s=%.3f%s', $stats_names[$i],$avg,"\n";
}
my $total_bid=$avg_all_files[0];
my $total_mono=$avg_all_files[1];
my $total_none=$avg_all_files[2];
my $AVOGADRO = 6.02e23;
my $VOLUME = 2.7e-23;
#my $kd = ($total_none*4)/(($total_mono+$total_bid)*$AVOGADRO*$VOLUME);
printf '%s=%.3f%s', "kd",$kd,"\n";


#subroutine gather : 
# returns percentage of bidentate, monodentate, no hydrogen bonds, and Dipeptide-Inositol-Inositol-Dipeptide, multiple inositol bound
# outputs the ramachandran map files for each of the bidentate, monodentate, none binding populations, and avg inositols bound.
# *.class file fields
# 0)NONE 1)CO0CO2_12 2)CO0NH3_12  3)NH1CO2_12  4)NH1NH3_12  5)CO0CO2_13  6)CO0NH3_13 7)NH1CO2_13 
#        8)NH1NH3_13 9)CO0CO2_14 10)CO0NH3_14 11)NH1CO2_14 12)NH1NH3_14 13)MONO    14)isDIID 
#	15)Phi      16)Psi

sub gather{
	open(FILE, "$filename");
	my ($num_bid,$num_mono,$num_none,$num_diid,$multi_inos,$sum_num_inos_bound,$mono_and_bid)=0;
	#my @bid_states;
	my $total=0;	# total number of snapshots (lines ) read
	
	while($line=<FILE>){
		@fields=split(/\s+/,$line);
		#count bidentate states
		my $sum1_12=0;
		for(my $i=1; $i<=12; $i++){
			$sum1_12+=$fields[$i];
		}
		if($sum1_12 > 0){
			#snapshot is a bidentate state
			#count bidentate
			$num_bid++;
			
			#output to one of 12 rama map files
			if($fields[1]){
				print BID1 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[2]){
				print BID2 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[3]){
				print BID3 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[4]){
				print BID4 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[5]){
				print BID5 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[6]){
				print BID6 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[7]){
				print BID7 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[8]){
				print BID8 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[9]){
				print BID9 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[10]){
				print BID10 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[11]){
				print BID11 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
			if($fields[12]){
				print BID12 "$filename:$total $fields[15] $fields[16] $fields[17]\n";
			}
		}	
		$total++;
	}
}
