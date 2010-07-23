#!/usr/bin/perl
if(@ARGV < 1){
	print "usage plot_sep_rama.pl <output-file-name>\n";
	exit;
}
@ramafiles=<*.dat>;
$output = $ARGV[0];
open(PLOT, ">${output}.p");
print PLOT "set size square\n";
print PLOT "set term post color enhanced 10\n";
print PLOT "set output \"${output}.ps\"\n";
print PLOT "set xrange [-180:180]\n";
print PLOT "set yrange [-180:180]\n";
print PLOT "set key bottom\n";
print PLOT "set multiplot\n";
print PLOT "set size 0.33,0.5\n";
$x=0;
$y=0;
foreach $file (@ramafiles){
	$option = "w dots";
	$size = -s $file;
	print $file, $size,"\n";
	if($size !=0){
		if($size < 5000){
			$option= "w p"; 
		}
		print PLOT "set origin $x, $y\n";
		print PLOT "plot \"$file\" $option\n";
	
		if($x==0.66){
			$x=0;
			$y+=0.5;
			if($y == 1){
				$y=0;
				print PLOT "unset multiplot\n";
				print PLOT "set multiplot\n";
			}
		} else {
			$x+=0.33;
		}
	}
}

print PLOT "unset multiplot\n";
