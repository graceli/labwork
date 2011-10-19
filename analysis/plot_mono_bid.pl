#!/usr/bin/perl

#script to generate a customized
#gnuplot script to produce a bidentate and monodentate plot in ps for the inositol alanine dipeptide system
if(@ARGV < 5) {
	print "usage: <x-scale> <y-scale> <simulation-length> <mono-file> <bid-file>\n";
	exit;
}

$xscal=$ARGV[0];
$yscal=$ARGV[1];
$simlen=$ARGV[2];
$monofile=$ARGV[3];
$bidfile=$ARGV[4];

$plot_fname="plot_ramas_inos_dipep.p";
open(PLOT, ">plot_ramas_inos_dipep.p");
print PLOT "set xlabel \"phi\"";
print PLOT "\n";
print PLOT "set ylabel \"psi\"";
print PLOT "\n";
print PLOT "set size square $xscal,$yscal";
print PLOT "\n";
print PLOT "set term post 10";
print PLOT "\n";
print PLOT "set out \"${monofile}.ps\"";
print PLOT "\n";
print PLOT "set title \"monodentate rama ($simlen ns)\"";
print PLOT "\n";
print PLOT "plot \"${monofile}\"";                  
print PLOT "\n";
print PLOT "set title \"bidentate rama ($simlen ns)\"";
print PLOT "\n";
print PLOT "set out \"${bidfile}.ps\"";              
print PLOT "\n";
print PLOT "plot [-200:200] [-200:200] \"${bidfile}\""; 

print PLOT "\n";
print PLOT "set out \"water_dipep_rama20ns.xvg.dat.ps\"";
print PLOT "\n";
print PLOT "set title \"alanine dipeptide rama ($simlen ns)\"";
print PLOT "\n";
print PLOT "plot \"water_dipep_rama20ns.xvg.dat\"";
print PLOT "\n";

#### multiplot part of script######
print PLOT "\n";
print PLOT "set out \"multi_mono_bid_rama.ps\"";
print PLOT "\n";
print PLOT "set multiplot";
print PLOT "\n";
print PLOT "set origin 0,0";
print PLOT "\n";
print PLOT "set title \"monodentate rama ($simlen ns)\"";
print PLOT "\n";
print PLOT "plot \"${monofile}\"";
print PLOT "\n";
print PLOT "set origin $xscal,0";
print PLOT "\n";
print PLOT "set title \"bidentate rama ($simlen ns)\"";
print PLOT "\n";
print PLOT "plot [-200:200] [-200:200] \"${bidfile}\"";
print PLOT "\n";
print PLOT "set title \"dipeptide in water ($simlen ns)\"";
print PLOT "\n";
print PLOT "set origin 0,$yscal";
print PLOT "\n";
print PLOT "plot [-200:200] [-200:200] \"water_dipep_rama20ns.xvg.dat\"";
print PLOT "\n";
print PLOT "unset multiplot\n";
system("gnuplot $plot_fname");

