#/usr/bin/perl
my $molname;

while($line=<STDIN>){
	my $title;
	if($line=~/[/ || $line =~/;/){
		#skip comment or title lines
		$title=$line;
		next;
	} else {
		# if line is blank ... I don't know how to match blank lines
		if($line =~ /[^0-9A-Za-z]/) {
			print "\n";
			next;
		}
		#otherwise
		if($title eq "[ atoms ]"){
			@f=split(/\s+/,$line);
			$atomNum = $f[0];
			$atomName = $f[1];
			$resNum = $f[2];
			$resName = $f[3];
			#4 is atomname again
			#5 is atom num again
			$charge = $f[6];
			$mass = $f[7];
			print "$atomNum  $atomName     $resNum  $resName  $atomName  $atomNum  $charge  $mass\n";
		}
	}
}


print "[ moleculetype ]\n";
print "; Name            nrexcl\n";
print "  $molname        3\n";

print "[ atoms ]\n";
print ";   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB\n";


