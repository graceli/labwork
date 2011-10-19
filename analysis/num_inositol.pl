#!/usr/bin/perl

# prints the number of waters in a gro file

if(@ARGV < 1){
	print "usage: num_waters <file>\n";
	exit;
}
$file=$ARGV[0];

chomp($num_waters=`cat $file | grep INS | wc -l`);

print $num_waters/24,"\n";

