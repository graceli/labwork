#!/usr/bin/perl
#
if(@ARGV != 1){
	print "Usage: perl grace2gnuplot.pl <grace-plot-name>\n";
	die;
}
my $graceplot = $ARGV[0];

system("sed \'s/@/#/g\' $graceplot > ${graceplot}.dat");

open(OUT, ">${graceplot}.p");

print OUT "plot \'${graceplot}.dat\' u 1:2\n";
