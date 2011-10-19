#!/usr/bin/perl
use strict;

if(@ARGV < 2){
	print "usage: bring_down.pl <dir_from> <cluster>\n";
	exit;
}

my $dir_from = $ARGV[0];
my $cluster = $ARGV[1];

print "Would you like to use dry-run?\n";
chomp(my $ans=<STDIN>);
my $use_dry_run="";
if($ans eq "y"){
	$use_dry_run = "--dry-run";
}

print "Would you like to use remove-sent-files?\n";
chomp(my $ans2=<STDIN>);
my $use_remove="";
if($ans2 eq "y"){
	$use_remove = "--remove-sent-files";
}

print "Would you like to get success dirs?\n";
chomp(my $ans3=<STDIN>);
my $get_success=$ans3;

#location of rsync on CCB 
my $use_rsync="";
if($cluster =~ /ccb/){
	$use_rsync="--rsh=ssh --rsync-path=/tools/rsync/2.6.9/bin/rsync";
}

system("rsync $use_rsync -av --progress $use_dry_run $use_remove ${cluster}:$dir_from/*/*/*.xtc .");

#print("rsync $use_ccb_rsync -av --progress $use_dry_run $use_remove ${cluster}:$dir_from/DATA/*/*.xtc .");

if($get_success eq "y"){
	system("rsync $use_rsync -av --progress $use_dry_run $use_remove ${cluster}:$dir_from/*_success/*desorted.xtc .");
}






