#!/usr/bin/env python

import glob
import numpy
from optparse import OptionParser

def error(data):
	pass

def aggregate(data):
	print data
	nrow,ncol = data.shape
	average=[]
	std=[]
	for i in range(0,ncol):
		average.append(numpy.sum(data[:,i],axis=0))
		std.append(numpy.std(data[:,i]))

	return average,std

if __name__ == "__main__":
	usage =  """ %prog [options] [shell expression] """
	parser = OptionParser(usage)
	#parser.add_option
	(options,args) = parser.parse_args()
	
	if len(args)<1:
		parser.error("Missing input files")

	files = args[0]	
	fileslist = glob.glob(files);


	data = numpy.array([])
	count = 0
	total_sampling = 3924038
	histogram = numpy.array([0]*7)

	for file in fileslist:
		A = numpy.genfromtxt(file)
		summedA = A.sum(axis=0)		
		histogram = numpy.add(histogram, summedA)

		count+=1
		if count == 100:
			#add A to a column in a numpy array called data
			print "added a row"
			print data.shape
			print histogram.shape
			data = numpy.append(data,(1.0/total_sampling)*histogram)
			print data
			histogram = numpy.array(numpy.arange(0,7))
			count = 0

	data = numpy.append(data,(1.0/total_sampling)*histogram)
	data.shape = (data.size/7,7)
	a,s = aggregate(data)

	# this is hardcoded

	print "average = ", a
	print "std = ", s
	
