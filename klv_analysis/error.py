#!/usr/bin/env python

import glob
import numpy
from optparse import OptionParser

def error(data):
	pass

# computes a list of averages and std
def aggregate(data, type="average"):
	print data
	nrow,ncol = data.shape
	average=[]
	std=[]
	for i in range(0,ncol):
		average.append(numpy.average(data[:,i],axis=0))
		std.append(numpy.std(data[:,i]))

	return numpy.array(average),numpy.array(std)

if __name__ == "__main__":

	# Implements block averaging over data that is stored as a row vector
	# Currently, this script reads in a list of files
	# Each file containing 2D matrix data with data from each frame
	# in a row

	usage =  """ %prog [options] [shell expression] """
	parser = OptionParser(usage)
	parser.add_option("-o", dest="outfilename", help="write to file", metavar="FILE")

	(options,args) = parser.parse_args()
	
	if len(args)<1:
		parser.error("Missing input files")

	files = args[0]	
	fileslist = glob.glob(files);


	data = numpy.array([])
	total_sampling = 3924038
	histogram = numpy.array([0]*7)

	block_rows=0
	count = 0
	for file in fileslist:
		A = numpy.genfromtxt(file)
		nrows, ncols = A.shape
		summedA = A.sum(axis=0)		
		histogram = numpy.add(histogram, summedA)

		count+=1
		block_rows += nrows
		if count == 1:
			#add A to a column in a numpy array called data
#			print "added a row"
#			print data.shape
#			print histogram.shape
#			print block_rows

			scale = 1.0/block_rows

			# don't apply scaling
			data = numpy.append(data,(1.0/100200/4)*histogram)
#			print data
			histogram = numpy.array(numpy.array([0]*7))
			count = 0
			block_rows = 0


	print "last histogram", histogram

#	data = numpy.append(data,histogram)
	data.shape = (data.size/7,7)
	a,s = aggregate(data)

	# this is hardcoded

	print "data"
	print data

	print "average = ", a
	print "std = ", s
	
	numpy.savetxt(options.outfilename+'.txt', numpy.transpose([a,s]), fmt="%0.3f %0.3f")	

