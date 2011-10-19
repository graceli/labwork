#!/usr/bin/perl
use strict;

if(@ARGV < 1){
	print "perl bind_pattern.pl <output-file>\n";
	exit;
}

open(FILE, "$ARGV[0]");
my @lines=<FILE>;  #read the output file lines into an array

#variables tracking 2 HB 2 groups binding
my $COCO=0;
my $CONH=0;
my $NHNH=0;
my $others=0;
my $mono=0;
#variables tracking 2HB single group binding
my $single_bid=0;
my $single_COCO=0;
my $single_CONH=0;
my $single_NHNH=0;
my $single_others=0;

for(my $i=0; $i<@lines; $i++){
	my @linefields = split(/\s+/,$lines[$i]);
	
	#split by the pipe to get info on each inositol
	#ie each $inosfield[$i] represents an inositol
	my @inosfields = split(/\|/,$linefields[6]);
	for(my $i=0; $i < @inosfields; $i++){
		#if inositol is bound to something count
		count_bidentate($inosfields[$i]);
		count_monodentate($inosfields[$i]);
		count_single_bid($inosfields[$i]);
	}
}

print "### 1 HB ###\n";
print "mono = $mono\n";
print "### 2 HB ###\n";
print "COCO = $COCO\n";
print "CONH = $CONH\n";
print "NHNH = $NHNH\n";
print "others = $others\n";
print "single-bid = $single_bid\n";
print "COCO = $single_COCO\n";
print "CONH = $single_CONH\n";
print "NHNH = $single_NHNH\n";
print "others = $single_others\n";

sub count_monodentate {
	my $inosfields = $_[0];
	if($inosfields){
		#get info on each groups of inositol
		my @ohfields1=split(/,/,$inosfields);
		#if @ohfields1 is 1, this means only 1 OH group is bound to something
		if(@ohfields1 == 1){
			#check if this field is empty -- not sure if it's really needed
			if($ohfields1[0]){
				#check that the OH group is only bound to 1 of something
				#there are cases that it is bound to more than 2 of something
				my @gfields = split(/;/,$ohfields1[0]);
				if(@gfields == 1){
					#print "$inosfields @gfields\n";
					$mono++;
				}
			}
		}
	}
}

#counts the number inositols forming 2 HB using 1 OH group
sub count_single_bid {
	my $inosfields = $_[0];
	if($inosfields){
		#check that it's only 1 bound group
		my @ohfields1=split(/,/,$inosfields);
		if(@ohfields1 == 1){
			if($ohfields1[0]){
				my @gfields = split(/;/,$ohfields1[0]);
				#make sure this group is bound to exactly 2 of something
				if(@gfields == 2){
					#print "$inosfields @gfields\n";
					$single_bid++;
					#since it's bound to 2 groups, classify it's specific
					#type
					my $a = $gfields[0];
					my $b = $gfields[1];
					if($a=~/CO/ && $b=~/CO/){
						$single_COCO++;
					}elsif($a=~/NH/ && $b=~/CO/ || $a=~/CO/ && $b=~/NH/){
						$single_CONH++;
					}elsif($a=~/NH/ && $b=~/NH/){
						$single_NHNH++;
					}else{
						$single_others++;
					}
				}
			}
		}
	}
}

#counts the number of inositols forming 2 HB with 2 different OH groups
#breaks down to COCO, CONH, NHNH
#this function implicitly assumes that a single inositol cannot form HBs using
#3 OH groups at once with a "flat" sheet structure. But once, the sheet starts deforming 
#significantly, this will not be true
sub count_bidentate {
	my $inosfields = $_[0];
	if($inosfields){
		my @ohfields1=split(/,/,$inosfields);
		if($ohfields1[0] && $ohfields1[1]){
			#print "@ohfields0\n";
			#split by ; (multiple bound groups per OH group)

			my @gfields0 = split(/;/,$ohfields1[0]);
			my @gfields1 = split(/;/,$ohfields1[1]);

			#print "$gfields0[0] $gfields1[0]\n";
			if(@gfields0 == 1 && @gfields1 == 1) {
				my $a = $gfields0[0];
				my $b = $gfields1[0];
				if($a=~/CO/ && $b=~/CO/){
					$COCO++;
				}elsif($a=~/NH/ && $b=~/CO/ || $a=~/CO/ && $b=~/NH/){
					$CONH++;
				}elsif($a=~/NH/ && $b=~/NH/){
					$NHNH++;
				}else{
					$others++;
				}
			}
		}
	}
}

