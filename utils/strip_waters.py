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

# figure out how to do this with subprocess
os.system("echo -e '!\"SOL\"\nq' | make_ndx -f %(tprfile)s -o nosol.ndx" % vars())
os.system("echo -e '!SOL' | trjconv -f %(tprfile)s -s %(tprfile)s -o nosol.gro -n nosol.ndx" % vars())

index=0
for xtcfile in xtcList:
	output=str(index)+"_nosol"
	if not os.path.exists(output):
		selection = '!SOL'
		command = "trjconv -f %s -s %s -n %s -o %s" % (xtcfile, tprfile, options.outputFilename, output)
		returnCode = subprocess.Popen(shlex.split(command),
		stdin=subprocess.PIPE).communicate("%s" % (selection))
		index += 1

