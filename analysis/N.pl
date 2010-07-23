#/usr/bin/perl

#given the concentration, and volume, determine teh number of molecules needed to make that
#concentration in solution

if(@ARGV < 2){
	print "usage: N.pl <conc> <volume nm cubed>\n";
	exit;
}

$concentration=$ARGV[0];
$volume=$ARGV[1];


$avogadro=6.02e23;
$nm2liter=1e-24;

$numMolec=$concentration*$avogadro*$nm2liter*$volume;
print $numMolec, "\n";


