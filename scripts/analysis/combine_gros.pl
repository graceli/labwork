#!/usr/bin/perl

# Last updated: Oct 2008 by Grace 
# this script takes two gro files and combines them into one
# warning:  the final combined gro file may not have the correct bounding box
# and will require the removal of the box dimensions from the final gro file and re-running editconf -f ... -c
# which will then put the correct bounding box

if(@ARGV < 1){
	print "usage: combine_gros.pl <gro1 file> <gro2 file> <output filename>\n";
	exit;
}
$gro1=$ARGV[0];
$gro2=$ARGV[1];
$output=$ARGV[2];

chomp($gro1b=`basename $gro1 .gro`);
chomp($gro2b=`basename $gro2 .gro`);

# convert the files into pdb
system("editconf -f $gro1 -o ${gro1b}.pdb");
system("editconf -f $gro2 -o ${gro2b}.pdb");

system("sed -e \'s/TER\$//\' -e \'s/ENDMDL\$//\' ${gro1b}.pdb > combo.pdb");

# what does this line do ???
system("sed -e :a -e \'/^\\n\*\$\/{\$d;N;\};\/\\n\$\/ba\' combo.pdb > combo2.pdb");

# combine pdbs (suppresse file name printing from grep)
system("grep -h ATOM ${gro2b}.pdb >> combo2.pdb");
system("echo TER >> combo2.pdb");
system("echo ENDMDL >> combo2.pdb");

# convert back into gro files
system("editconf -f combo2.pdb -o $output");
