#!/usr/bin/perl

#this script is pretty much deprecated since genconf can make a box of peptides from a trajectory
# originally called translate.pl

# select 4 random chains from the folder originals/
# and add them to the corners of a user specified box, such that 
# peptides are pairwise equidistant
# The big box size is 30 angstroms cubed
# The small box is determined to be 15 angstroms cubed
# reusing the code from add_inositol.pl

if(@ARGV < 3) {
	print "perl translate.pl <peptide-name> <number of atoms> <comment>\n";
	exit;
}



$peptide_name = $ARGV[0];
$number_atoms = $ARGV[1];
$comment = $ARGV[2];

if($peptide_name eq "klvffae") {
	$struct_dir="\/work\/grace\/traj\/structs\/klvffae";
}

if($peptide_name eq "f8") {
	$struct_dir="\/work\/grace\/traj\/structs\/f8";
}

if($peptide_name eq "q8") {
	$struct_dir="\/work\/grace\/traj\/structs\/q8";
}

@smallbox=(4,4,4);
@sbo=(0,0,0);

$a = int(rand(500));
$b = int(rand(500));
$c = int(rand(500));
$d = int(rand(500));

print STDERR "$a $b $c $d\n";

$total_atoms=4*$number_atoms;

print "$peptide_name peptides used ($a,$b, $c, $d); $comments\n";
print "  $total_atoms\n";

translmmp($a, $number_atoms);
translpmm($b, $number_atoms);
translmpm($c, $number_atoms);
translppp($d, $number_atoms);


#sub routines

sub translmmp {
        my $a=$_[0];
        my $total_atoms = $_[1];

        print STDERR "$total_atoms\n";

	open(FILE, "$struct_dir/${peptide_name}_${a}.gro") || die "not opened";

	@file = <FILE>;
        $last_line = $total_atoms + 1;
	for $i (2..${last_line}){
		#remove blanks at the beginning of line
		$file[$i]=~s/^\s+//;
		@f = split(/\s+/,$file[$i]);
		#translate by 0.75 nm
		$f[3]= $f[3] - $smallbox[0]/2 - $sb0[0];
		$f[4]= $f[4] - $smallbox[1]/2 - $sb0[1];
		$f[5]= $f[5] + $smallbox[2]/2 + $sb0[2];
		$f[0]=~m/(\d+)(.+)/;
		#print $1,$2;
		printf '%5d%-5s%5s%5d%8.3f%8.3f%8.3f',$1,$2,@f[1..5];
		#printf '%5d%-5s',$1,$2;
		print "\n";
	}
	close(FILE);
}

sub translmpm{
        my $a=$_[0];
        my $total_atoms = $_[1];

        print STDERR "$total_atoms\n";

	open(FILE, "$struct_dir/${peptide_name}_${a}.gro");
	@file = <FILE>;
        $last_line = $total_atoms + 1;

	for $i (2..${last_line}){
		#remove blanks at the beginning of line
		$file[$i]=~s/^\s+//;
		@f = split(/\s+/,$file[$i]);
		#translate by 0.75 nm
		$f[3]= $f[3] - $smallbox[0]/2 - $sb0[0];
		$f[4]= $f[4] + $smallbox[1]/2 + $sb0[1];
		$f[5]= $f[5] - $smallbox[2]/2 - $sb0[2];
		$f[0]=~m/(\d+)(.+)/;
		#print $1,$2;
		printf '%5d%-5s%5s%5d%8.3f%8.3f%8.3f',$1,$2,@f[1..5];
		#printf '%5d%-5s',$1,$2;
		print "\n";
	}
	close(FILE);
}

sub translpmm{
        my $a=$_[0];
        my $total_atoms = $_[1];

        print STDERR "$total_atoms\n";

	open(FILE, "$struct_dir/${peptide_name}_${a}.gro");
	@file = <FILE>;
        $last_line = $total_atoms + 1;
	
        for $i (2..${last_line}){
		#remove blanks at the beginning of line
		$file[$i]=~s/^\s+//;
		@f = split(/\s+/,$file[$i]);
		#translate by 0.75 nm
		$f[3]= $f[3] + $smallbox[0]/2 + $sb0[0];
		$f[4]= $f[4] - $smallbox[1]/2 - $sb0[1];
		$f[5]= $f[5] - $smallbox[2]/2 - $sb0[2];
		$f[0]=~m/(\d+)(.+)/;
		#print $1,$2;
		printf '%5d%-5s%5s%5d%8.3f%8.3f%8.3f',$1,$2,@f[1..5];
		#printf '%5d%-5s',$1,$2;
		print "\n";
	}
	close(FILE);
}

sub translppp{
        my $a=$_[0];
        my $total_atoms = $_[1];
        
        print STDERR "$total_atoms\n";

        open(FILE, "$struct_dir/${peptide_name}_${a}.gro");
	@file = <FILE>;
        $last_line = $total_atoms+1;    
    
	for $i (2..${last_line}){
		#remove blanks at the beginning of line
		$file[$i]=~s/^\s+//;
		@f = split(/\s+/,$file[$i]);
		#translate by 0.75 nm
		$f[3]= $f[3] + $smallbox[0]/2 + $sb0[0];
		$f[4]= $f[4] + $smallbox[1]/2 + $sb0[1];
		$f[5]= $f[5] + $smallbox[2]/2 + $sb0[2];
		$f[0]=~m/(\d+)(.+)/;
		#print $1,$2;
		printf '%5d%-5s%5s%5d%8.3f%8.3f%8.3f',$1,$2,@f[1..5];
		#printf '%5d%-5s',$1,$2;
		print "\n";
	}
	close(FILE);
}



