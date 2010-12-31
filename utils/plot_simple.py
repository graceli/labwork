#!/usr/bin/env python

import pylab
import os
from optparse import OptionParser

usage = "%prog [options] [flat file]"

parser = OptionParser(usage)
options,args = parser.parse_args()

if len(args) < 1:
	parser.error("Missing input file")

txtfile = args[0]

if not os.path.exists(txtfile):
	print "Error:", file, "does not exist"
	sys.exit(1)


dataArray = pylab.genfromtxt(txtfile)
nrows, ncols = dataArray.shape

x = dataArray[:,0]
y = dataArray[:,1]

pylab.plot(x,y)
pylab.axis('tight')
pylab.grid(True)
pylab.savefig('fig.png')
pylab.show()


