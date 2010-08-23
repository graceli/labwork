#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import shlex
import glob
import sys

from optparse import OptionParser

usage = """%prog [options] <tprfile> <xtcfile> """

parser = OptionParser(usage)

parser.add_option("-o", dest="outputFilename", default="nosol",
help="specify the output trajectory file name")

(options, args) = parser.parse_args()

if len(args)<2:
	parser.error("Missing input files")

tprfile = args[0]
files = args[1]

if not os.path.exists(tprfile):
	print "Error: the tpr file was not found"
	sys.exit(0)

if len(files) < 1:
	print "Error: the xtc file(s) was not specified"
	sys.exit(0)

print files
xtcList = glob.glob(files)
print xtcList

#create a no nosolvent index file
# I can't get make_ndx to read the \n properly here
#selection = '!"SOL"'
#command = "make_ndx -f %s -o %s" % (tprfile, options.outputFilename)
#returnCode = subprocess.Popen(shlex.split(command), stdin=subprocess.PIPE).communicate("%s%s%s" % (selection,"\n", 'q'))
#print returnCode
#sys.exit(0)

index=0
for xtcfile in xtcList:
	selection = '!SOL'
	command = "trjconv -f %s -s %s -n %s -o %s" % (xtcfile, tprfile, options.outputFilename, str(index)+"_nosol")
	returnCode = subprocess.Popen(shlex.split(command),
	stdin=subprocess.PIPE).communicate("%s" % (selection))
	index += 1


