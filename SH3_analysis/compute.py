import tables
import numpy

def main():
	"""docstring for main"""
	h5file = tables.openFile('DR_analysis.h5')
	
	
	# data_640K = h5file.read(where='0.792 > w >= 0.786')['w']
	data_300K = h5file.read(where='1.677 > w >= 1.689')['w']
	
	average = numpy.average(data_300K)
	
	print "<rg> = ", data_300K
	
	
	# [Mar/01/2011 00:37:45]        1.677      500000           3      24297.379        300.000           0.935
	# [Mar/01/2011 00:37:45]        1.689      500000           3      22693.467        298.000           0.923