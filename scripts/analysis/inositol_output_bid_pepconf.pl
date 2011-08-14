#!/usr/bin/perl
use warnings;
use strict;
use POSIX qw(ceil floor);

# outputs a files: bidentate, monodentate and none. Each line of the file containing a number 0-2 indicating the binding state of inositol, and the (phi, psi) angles of the alanine dipeptide

my @class=<*.class>;

#process each *.class file
foreach my $cfile (@class){
	gather($cfile);	
}	

sub gather{
	my $filename = $_[0];
	open(FILE,"$filename");
	my $total=0;	# total number of snapshots (lines ) read
	my @bid;
	my $line;
	my @fields;	
	while($line=<FILE>){
		@fields=split(/\s+/,$line);

		#count bidentate states
		my $sum1_12=0;
		my $phi = $fields[15];
		my $psi = $fields[16];

		for(my $i=1; $i<=12; $i++){
			$sum1_12+=$fields[$i];
		}
		if($sum1_12 > 0){
			#snapshot is a bidentate state
			print "2 $phi $psi\n";
		}elsif ($fields[13] > 0) {
			print "1 $phi $psi\n";
		}elsif ($fields[0] == 4){
			print "0 $phi $psi\n";
		}
		$total++;
	}
}

sub DEBUG{
	print STDERR $_[0],"\n";
}
