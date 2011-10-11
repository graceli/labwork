#!/usr/bin/env python
import numpy
import sys
import pylab
#import config
import subprocess

def reorder_and_normalize(map_filename, data):
	""" reorder the nonpolar matrix columns using the mapping in a flat file called table.dat 
		mapping is a set of (key, value) pairs where key = residue and value = number of atoms
		this file is only used for inositol-nonpolar residue binding and is not sorted by key
	"""
	import csv
	import re
	
	#load in residues
	r=csv.reader(open(map_filename), delimiter=' ')
	residue_map=[]
	residue_map_dict = {}
	for line in r:
		residue_map.append(line[0])
		residue_map_dict[line[0]]=int(line[1])

	nonpolar_dict={}
	# print data.shape
	# build dictionary 
	for i in range(0, len(residue_map)):
		# print i
		nonpolar_dict[residue_map[i]] = data[:,i]

	# print "non-sorted data", nonpolar_dict
	
	data_in_order = []
	for key in sorted(nonpolar_dict.keys(),key=lambda k: int(re.search('[0-9]+', k).group(0))):
		print "storing key", key
		#return data in the order of residue number and normalize by the number of atoms 
		data_in_order.append(nonpolar_dict[key]/int(residue_map_dict[key]))

	return numpy.transpose(numpy.array(data_in_order))

def smooth(data, window_size, time_present=True, timestep=1):
	if len(data.shape) > 1:
		nrows,ncols = data.shape
		assert nrows > ncols, "Number of cols is more than rows in the matrix!"
	else:
		nrows = data.shape[0]
		ncols = 1

	window_actual = window_size/timestep
	total_length = nrows - window_actual + 1
	smoothed = []
	for w in range(0, total_length):
		if time_present == False:
			values = [ w*timestep + 0.5*window_size ]
		else:
			values = []

		if ncols == 1:
			values.append(numpy.average(data[w:w+window_actual + 1]))
		else:
			for i in range(0, ncols):
				values.append(numpy.average(data[w:w+window_actual + 1, i]))
	
		smoothed.append(values)

	if smoothed == []:
		print "WARNING: smoothed is empty"
		return smoothed
	else:
		return numpy.array(numpy.vstack(smoothed))

def columnAverage(data,colnum):
	total_num_systems = len(data)
	stats = []
	ndatapoints = len(data[0][:,1])
	print "there are",ndatapoints, "datapoints" 
	for time in range(0, ndatapoints):
		values = []
		for i in range(0, total_num_systems):	
			if time < len(data[i][:,colnum]):
				values.append(data[i][time,colnum])

		avg = numpy.average(values)
		std = numpy.std(values)
		stats.append([avg, std])	

	print "computed", len(stats), "number of points"
	return numpy.array(stats)
		
def savetxt(filename, header, data, compress='.gz', fmt='%0.2f'):
	f=open(filename, 'w')
	f.write(header)
	f.close()
	numpy.savetxt(filename, data, fmt=fmt)
	# retcode = subprocess.call(["gzip", filename])
	# f.close()

def summary_statistics(data, sum_across="columns"):
	# hack -- data matrix format is specific to KLVFFAE
	# data -- list of time series and becomes rows in a matrix when vstacked
	# transpose is done to make time series run up and down, which is more intuitive to me
	print data
	# data = numpy.transpose(numpy.vstack(data))
	nrows, ncols = data.shape
	print data

	assert nrows > ncols
	axis=0
	if sum_across == "columns":
		axis=1

	return [numpy.average(data, axis=axis), numpy.std(data, axis=axis)]

def array_list_to_matrix(arraylist):
	mat = numpy.transpose(numpy.vstack(numpy.array(arraylist)))
	nrows, ncols = mat.shape
	
	assert nrows > ncols, "Expected matrix to have more rows than columns"
	
	return mat

def main():
	"""docstring for main"""
	file = sys.argv[1]
	data = numpy.genfromtxt(file)
	print data

	data_smoothed = smooth(data, 500)
	numpy.savetxt(file+"_smoothed.txt", data_smoothed)

	#print data_smoothed
	
if __name__ == '__main__':
	main()
