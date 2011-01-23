#!/usr/bin/env python
import numpy
import sys

def smooth(data, window_size):
	nrows,ncols = data.shape
	total_length = nrows - window_size + 1
	smoothed = []
	print smoothed
	for w in range(0, total_length):
		values = [ w + 0.5*window_size ]
		for i in range(0, ncols):
			values.append(numpy.average(data[w:w+window_size+1, i]))
		smoothed.append(values)
	return numpy.array(smoothed)

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