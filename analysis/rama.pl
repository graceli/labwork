#!/usr/bin/perl
use strict;

#read 4 hbondnum files to determine whether a snapshot is bidentate or monodentate
#if (snapshot bidentate)
#	read rama file for dihedral angles
#	output to bidentate rama plot
#else if (snapshot monodentate)
#	read rama file for dihedral angles
#	output to bidentate rama plot
my($file4,$file5,$file6,$file7);
my(@inos4,@inos5,@inos6,@inos7);

if(@ARGV<5){
	print "usage: <hbnum1> <hbnum2> <hbnum3> <hbnum4> <rama-file-total>\n";
	exit;
}
$file4 = $ARGV[0];
$file5 = $ARGV[1];
$file6 = $ARGV[2];
$file7 = $ARGV[3];
my $rama = $ARGV[4];

@inos4 = read_xvg_file($file4);
@inos5 = read_xvg_file($file5);
@inos6 = read_xvg_file($file6);
@inos7 = read_xvg_file($file7);
my @rama = read_xvg_file($rama);

if (scalar(@rama) == scalar(@inos4)){
	print STDERR "size of rama equals size of inos4\n";
}


my $i;

open(BID, ">${rama}_bid_rama.dat");
open(MONO,">${rama}_mono_rama.dat");

for($i=0; $i<@inos4; $i++){
	my $hbnum4=gethbnum($inos4[$i]);
	my $hbnum5=gethbnum($inos5[$i]);
	my $hbnum6=gethbnum($inos6[$i]);
	my $hbnum7=gethbnum($inos7[$i]);


	if( $hbnum4 == 2 || $hbnum5 == 2 || $hbnum6 == 2 || $hbnum7 == 2){
		#bidentate
		print "bid\n";
		print BID $rama[$i];	
	}elsif ($hbnum4 == 1 || $hbnum5 == 1 || $hbnum6 == 1 || $hbnum7 == 1){
		print "mono\n";
		print MONO $rama[$i];
	}else{
		print "none\n";
	}
}

sub gethbnum{
	return unpack("x21 A1",$_[0]);
}

sub read_xvg_file{
	my $file=$_[0];
	my @contents;
	my $line;
	open(HB, $file);
	while($line=<HB>){
		if($line =~/#/ || $line =~ /@/) {
			next;
		}
		push(@contents,$line);
	}

	return @contents;
}

