#!/usr/bin/perl
#Grace Li October 8, 2007
#This script is a wrapper for qsub -l kind -b y 
#Eliminates the need to cd to current directory
#Currently, to use with Gromacs binaries, must send in
#full path to where the binaries are located

if(@ARGV < 1){
	print "usage: kind \"<command-to-execute>\"\n see kind -h for examples\n";
	exit;
}

if($ARGV[0] eq "-h"){
	print "Usage examples for kind:\n";
	print "single command:\n";
	print "\tkind \"gunzip foo.tar.gz\"\n";
        print "multiple commands:\n";
	print "\tkind \"gunzip foo.tar.gz; tar xvvf foo.tar; cd foo\"\n";
	print "Gromacs commands:\n";
	print "\tkind \"/full/path/to/grompp -f empty.mdp -c out.gro -p out.top -o out.tpr\"\n";
	exit;
}

$pwd=$ENV{'PWD'};
system("qsub -l kind -N kind -b y \"cd $pwd;$ARGV[0]\"");

