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

DEBUG("### $filename ###:");
DEBUG("### calculating error for bidentate ###");
my $state = "bid";
my @a;
@a = gather($filename);
error();

DEBUG("### calculating error for monodentate ###");
$state = "mono";
@a = gather_mono($filename);
error();

DEBUG("### calculating error for none ###");
$state = "none";
@a = gather_none($filename);
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
#returns an array of 0s and 1s corresponding to bidentate or not
sub gather{
	open(FILE,"$filename");
	my $total=0;	# total number of snapshots (lines ) read
	my @bid;
	my $line;
	my @fields;	
	while($line=<FILE>){
		@fields=split(/\s+/,$line);

		#count bidentate states
		my $sum1_12=0;
		for(my $i=1; $i<=12; $i++){
			$sum1_12+=$fields[$i];
		}
		if($sum1_12 > 0){
			#snapshot is a bidentate state
			#count bidentate
			push(@bid,1);
		}else{
			push(@bid,0);
		}
		$total++;
	}
	return @bid;
}

sub gather_mono{
	open(FILE,"$filename");
	my $total=0;	# total number of snapshots (lines ) read
	my @mono;	
	my $line;
	my @fields;
	while($line=<FILE>){
		@fields=split(/\s+/,$line);

		#count bidentate states
		if($fields[13] > 0){
			#snapshot is a monodentate state
			#count monodentate
			push(@mono,1);
		}else{
			push(@mono,0);
		}
		$total++;
	}
	return @mono;
}

sub gather_none{
	open(FILE,"$filename");
	my $total=0;	# total number of snapshots (lines ) read
	my @none;	
	my $line;
	my @fields;
	while($line=<FILE>){
		@fields=split(/\s+/,$line);
		#DEBUG("$fields[0]");
		if($fields[0] == 4){
			push(@none,1);
		}else{
			push(@none,0);
		}
		$total++;
	}
	return @none;
}
sub DEBUG{
	print STDERR $_[0],"\n";
}
