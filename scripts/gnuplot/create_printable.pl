#!/usr/bin/perl

# Grace Li
# Jan 2009
# this script takes a data file and generates a printable plot

use strict;

my $file = $ARGV[0];  #name of data file
my $plotname = $ARGV[1]; #name of plot file

open(OUT, ">${plotname}.p");
print OUT "set term post font 8\n";
print OUT "set output \"${plotname}.ps\"\n";
print OUT "set grid\n";
print OUT "set size 0.5, 0.5\n";
print OUT "plot \'$file\' u 1:(\$2-\$3):(\$2+\$3) with filledcurve lc rgb \"#cccccc\" notitle, \'\'  w l\n";
close(OUT);

system("'gnuplot' ${plotname}.p");

