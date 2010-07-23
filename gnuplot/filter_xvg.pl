#!/usr/bin/perl

@xvg=<rama/*.xvg>;

foreach $file (@xvg){
	open(FILE, "$file");
	open(OUT, ">${file}.fil");
	while($line=<FILE>){
		if($line=~/@/ || $line=~/#/){
			next;
		}
		chomp($line);
		@fields = split(/\s+/, $line);
		print OUT $fields[0]," ",$fields[1],"\n";
	}
	close(OUT);
	close(FILE);
}
