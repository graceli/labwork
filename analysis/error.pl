#!/usr/bin/perl
use warnings;
use POSIX qw(ceil floor);

open(IN, "bidentate_values");
open(OUT, ">s_plot2");
chomp(@a=<IN>);

my $block_size=1;
my $sum_blocks=0;

for(my $i=0; $i<@a;$i++){
	$sum_blocks+=$a[$i];
}
my $run_avg=$sum_blocks/@a;
print "run average = $run_avg\n";
my $run_var;

#for each block, of size $block_size
for(my $block_size=1; $block_size<=10; $block_size+=1){
	#for each block of size $block_size
	#outer loop sets the start of the next block
	#inner loop, loops over the entries of the block
	my @block_avg;
	my $block_var=0;
	my $num_blocks = floor(@a/$block_size);
	print "num_blocks = $num_blocks\n";
	for(my $i=0; $i<@a; $i+=$block_size){
		if(@a-$i < $block_size){
			last;
		}
		#compute block average
		my $block_sum=0;
		for(my $j=0; $j<$block_size; $j++){
			$block_sum+=$a[$i+$j];
		}
		my $block_avg = $block_sum/$block_size;
		$block_var += ($block_avg - $run_avg)**2;
		push(@block_avg,$block_avg);
	}
	$block_var=$block_var/$num_blocks;
	$s=$block_size*$block_var/0.40;		#run variance is 0.40
	#print "s=$s\n";
	my $actual_size = $block_size;
	my $s_actual = $s;
	my $inv_block_var = $block_var;
	print OUT "$block_size $s \n";
}
