#!/usr/bin/perl
my $ind = 1;
open(FIN,$ARGV[0]);
while (<FIN>)
{
	chomp;	
	printf ("%d %s\n",$ind,$_);
	$ind++;
}
close FIN;
