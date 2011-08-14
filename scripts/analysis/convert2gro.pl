#!/usr/bin/perl

print "usage: convert2gro.pl <isomer>\n";

my $num_threads=4;
my $residue="ala";
my $isomer=$ARGV[0];
my @t = (165,93,203,185,7);
my @seg = (5,10,10,10,10);
for(my $i=0; $i<5; $i++){
	for(my $s=0; $s<$seg[$i]; $s++){
		#convert each xtc into gros
		#my $gro = "echo 15 | trjconv -f ../xtc/a4${isomer}${t[$i]}_${s}.xtc -o ../gro/a4${isomer}${t[$i]}_${s}_part.gro -s ../tpr/a4${isomer}${t[$i]}_${s}.tpr -n /home/grace/work/inositol/index/ala_scy.ndx";  
		#produce rama files
		#my $rama = "g_rama -f ../xtc/a4${isomer}${t[$i]}_${s}.xtc -s ../tpr/a4${isomer}${t[$i]}_${s}.tpr -o rama/a4${isomer}${t[$i]}_${s}.xvg";
		#classification on each (gro,rama) pair
		my $class = "dipep_class ../gro/a4${isomer}${t[$i]}_${s}_part.gro rama/a4${isomer}${t[$i]}_${s}.xvg 2>/dev/null > a4${isomer}${t[$i]}_${s}.class &";

		#print "$gro\n";
		#print "$rama\n";
		print "$class\n";

		if($ARGV[1] eq "-exec"){
			#system($gro);
			#system($rama);
			system($class);
		}
	}
}
