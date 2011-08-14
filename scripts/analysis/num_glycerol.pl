#!/usr/bin/perl
if(@ARGV < 1){
	print "usage: num_octanes.pl <file>\n";
	exit;
}
chomp($num_lines=`cat $ARGV[0] | grep CER | wc -l`);
$num_octane=$num_lines/14;
print "$num_octane\n";
