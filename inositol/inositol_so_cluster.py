#!/usr/bin/python

import numpy

def average_cluster_columns(datafile, outputname):
	""" this function assumes that every odd column
		starting from 1 is time or is not used"""

	data = numpy.genfromtxt(datafile)
	sdata = data[:,1::2]
	print sdata
	averaged_columns = numpy.average(sdata, axis=1)
	print average_columns
	numpy.savetxt(outputname, numpy.transpose([data[:,0],averaged_columns]), fmt='%0.2f')


average_cluster_columns('scyllo_nclust_all.dat', 'scyllo.txt')
average_cluster_columns('chiro_nclust_all.dat', 'chiro.txt')
