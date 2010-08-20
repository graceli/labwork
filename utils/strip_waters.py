#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import shlex
from optparse import OptionParser

usage = """%prog [options] <tprfile> <xtcfile> """

 #usage = """
        #usage: %prog [options] <PSF file> <PDB file> <DCD file>
    #"""
    
    # initialize the parser with our custom usage string (above)
    #parser = optparse.OptionParser(usage)
    ## add all the command-line options
    #parser.add_option("-o", dest="h5_filename", default="analysis.h5",
#help="Output analysis file [default: %default]")    
    #parser.add_option("-s", dest="selection", default="name CA", help="Atom
#selection [default: %default]")    
    #parser.add_option("-t", dest="first_timestep", default=0, type="int",
#help="Starting timestep [default: %default]")    
    #parser.add_option("-f", dest="dcdtime", default=1, type="int", help="DCD
#output frequency [default: %default]")    
    #parser.add_option("-d", dest="dt", default=0.002, type="float",
#help="Integration Timestep [default: %default]")    
    
    ## parse the command line
    #(options, args) = parser.parse_args()
    
    ## check to make sure all arguments are there
    #if len(args) < 3:
        #parser.error("No input files specified")
    
    #psf_file = args[0]
    #pdb_file = args[1]
    #dcd_file = args[2]


parser = OptionParser(usage)

parser.add_option("-o", dest="outputFilename", default="nosol.xtc",
help="specify the output trajectory file name")

(options, args) = parser.parse_args()

if len(args)<2:
	parser.error("Missing input files")

tprfile = args[0]
xtcfile = args[1]

if not os.path.exists(tprfile):
	print "Error: the tpr file was not found"
	sys.exit(0)

if not os.path.exists(xtcfile):
	print "Error: the xtc file was not found"
	sys.exit(0)


#create a no nosolvent index file
selection = '!"SOL"'
command = "make_ndx -f %s -n %s" % (tprfile, outputIndex)
returnCode = subprocess.Popen(shlex.split(command),
stdin=subprocess.PIPE).communicate("%s" % (selection))

selection = '!SOL'
command = "trjconv -f %s -n %s -o %s" % (tprfile, xtcfile,
parser.outputFilename)
returnCode = subprocess.Popen(shlex.split(command),
stdin=subprocess.PIPE).communicate("%s" % (selection))