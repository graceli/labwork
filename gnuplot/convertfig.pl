#!/usr/bin/perl
if ($#ARGV<2)
{
  print "args input resolution output";
  exit;
}
system("gs -r$ARGV[1] -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=$ARGV[2] -dBATCH -dNOPAUSE $ARGV[0]");
