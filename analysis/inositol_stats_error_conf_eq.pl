#!/usr/bin/perl
use warnings;
use strict;
use POSIX qw(ceil floor);

if(@ARGV<4){
	print "Take a *.class file as input and produces the s vs block size plots for bidentate, mono , and none binding states\n";
	print "To determine the statistical inefficiency s, read off of the asymptote after plotting each s vs. block size\n";
	print "The error = sqrt(variance in the mean) = s*run_variance/length_of_run\n";
	print "usage: inositol_stats <filename.class> <start> <end> <step>\n";
	exit;
}
my $filename=$ARGV[0];

if(! -e $filename){
	print "$filename does not exist.\n";
	exit;
}
	
my $block_step=$ARGV[3];
my $start=$ARGV[1];
my $end=$ARGV[2];

my @R1_phi = (-180,-125); my @R1_psi = (100,180);
my @R2_phi = (-125, -25); my @R2_psi = ( 90,180);
my @R3_phi = (-180,-125); my @R3_psi = (-75,75);
my @R4_phi = (-120, -25); my @R4_psi = (-75,50);
my @R5_phi = (  50, 100); my @R5_psi = (-40,50);



DEBUG("### $filename ###:");
#DEBUG("### calculating error for bidentate and R1 ###");
my $state = "bid_R1";
my @a;
@a = gather(2,@R1_phi,@R1_psi);
error();

#DEBUG("### calculating error for monodentate and R1 ###");
$state = "mono_R1";
@a = gather(1,@R1_phi,@R1_psi);
error();

#DEBUG("### calculating error for none and R1 ###");
$state = "none_R1";
@a = gather(0,@R1_phi,@R1_psi);
error();

DEBUG("### calculating error for bidentate and R2 ###");
$state = "bid_R2";
@a = gather(2,@R2_phi,@R2_psi);
error();
DEBUG("### calculating error for monodentate and R2 ###");
$state = "mono_R2";
@a = gather(1,@R2_phi,@R2_psi);
error();
DEBUG("### calculating error for none and R2 ###");
$state = "none_R2";
@a = gather(0,@R2_phi,@R2_psi);
error();


DEBUG("### calculting error for bidentate and R3 ###");
$state = "bid_R3";
@a = gather(2,@R3_phi,@R3_psi);
error();
DEBUG("### calculating error for monodentate and R3 ###");
$state = "mono_R3";
@a = gather(1,@R3_phi,@R3_psi);
error();
DEBUG("### calculating error for none and R3 ###");
$state = "none_R3";
@a = gather(0,@R3_phi,@R3_psi);
error();

DEBUG("### calculating error for bidentate and R4 ###");
$state = "bid_R4";
@a = gather(2,@R4_phi,@R4_psi);
error();
DEBUG("### calculating error for monodentate and R4 ###");
$state = "mono_R4";
@a = gather(1,@R4_phi,@R4_psi);
error();
DEBUG("### calculating error for none and R4 ###");
$state = "none_R4";
@a = gather(0,@R4_phi,@R4_psi);
error();


DEBUG("### calculating error for bidentate and R5 ###");
$state = "bid_R5";
@a = gather(2,@R5_phi,@R5_psi);
error();
DEBUG("### calculating error for monodentate and R5 ###");
$state = "mono_R5";
@a = gather(1,@R5_phi,@R5_psi);
error();
DEBUG("### calculating error for none and R5 ###");
$state = "none_R5";
@a = gather(0,@R5_phi,@R5_psi);
error();


sub error{
	my ($sum_run,$sum_run_var,$s)=0;
	my $size = @a;
	DEBUG("total number of snapshots = $size");
	
	my $sum_blocks=0;
	open(OUT, ">s_plot2_${state}");
	
	#calculate the run average
	for(my $i=0; $i<@a;$i++){
		$sum_run+=$a[$i];
	}
	my $run_avg=$sum_run/@a;
	DEBUG("total = $sum_run");
	DEBUG("run average = $run_avg");
	
	#calculate the run variance
	# (1/size(a))*sum over all values of @a (a(i)-a_avg)^2
	for(my $i=0; $i<@a; $i++){
		$sum_run_var+=($a[$i]-${run_avg})**2
	}
	my $run_var=$sum_run_var/@a;
	DEBUG("run variance = $run_var");
	#DEBUG("block_start=$start  block_end=$end   block_step=$block_step");
	
	#outmost loop sets the $block_size
	#DEBUG("num blocks	block size");
	for(my $block_size=$start; $block_size<=$end; $block_size+=$block_step){	
		#dynamic array storing block avgs
		my @block_avg;
	
		#running sum of (block_avg-run_avg)^2
		my $block_var=0;
	
		#total number of blocks
		my $num_blocks = floor(@a/$block_size);
	
		#DEBUG("$num_blocks	$block_size");
		#This loop sets the start of a block of $block_size
		#essentially loops over each block
		for(my $i=0; $i<@a; $i+=$block_size){
			#if the number of elements left is smaller than a block size,skip out of the loop
			if(@a-$i < $block_size){
				#DEBUG("Not enough left to make a block...stopping");
				last;
			}
			
			#compute block average
			my $block_sum=0;
	
			#this loop, loops over each entry of the block
			for(my $j=0; $j<$block_size; $j++){
				$block_sum+=$a[$i+$j];
			}
			#compute the block mean (or average)
			my $block_mean = $block_sum/$block_size;
			#update the variance 
			$block_var += ($block_mean - $run_avg)**2;
			#store currently calculated block mean
			push(@block_avg,$block_mean);
		}
	
		#compute the variance in the block mean
		$block_var=$block_var/$num_blocks;
		#compute the statistical inefficiency for the simulation divided into blocks of $block_size
		$s=$block_size*$block_var/$run_var;
		my $sqrt_block_size=sqrt($block_size);
		print OUT "$block_size $s \n";
	}
}

#returns an array of 0s and 1s corresponding to bidentate AND PEP_CONF or not
# where PEP_CONF = R1, R2, ..., R5
sub gather{
	open(FILE,"$filename");
	my $total=0;	# total number of snapshots (lines ) read
	my @array;
	my $line;
	my @fields;	
	my $gather_state = $_[0];
	my @phi_range = ($_[1], $_[2]);
	my @psi_range = ($_[3], $_[4]);

	DEBUG("@phi_range @psi_range\n");

	while($line=<FILE>){
		@fields=split(/\s+/,$line);
		#count bidentate states
		my $state=$fields[0];
		my $phi = $fields[1];
		my $psi = $fields[2];
		#get binned phi, psi
		my $bin_phi = floor(($phi+180)/5)*5-180;
		my $bin_psi = floor(($psi+180)/5)*5-180;

		#count bidentate AND PEP_CONF
		if($state == $gather_state && ($phi >= $phi_range[0] && $phi < $phi_range[1] && $psi >= $psi_range[0] && $psi < $psi_range[1])) {
				push(@array,1);
		}else{
			push(@array,0);
		}
		
		$total++;
	}

	DEBUG("$total");

	return @array;
}

sub DEBUG{
	print STDERR $_[0],"\n";
}
