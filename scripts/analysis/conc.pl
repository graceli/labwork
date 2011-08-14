#/usr/bin/perl

#given the number of molecules determine the concentration in the given volume
# in M (mols per liter)

if(@ARGV < 2){
	print "usage: conc.pl <num mol> <volume nm cubed>\n";
	exit;
}

$numMolec=$ARGV[0];
$volume=$ARGV[1];
$avogadro=6.02e23;
$nm2liter=1e-24;
$concentration=$numMolec/($avogadro*$nm2liter*$volume);


print $concentration, "\n";
