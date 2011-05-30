import tables
import numpy

def main():
	"""docstring for main"""
	h5file = tables.openFile('DR_analysis.h5')
	print "reading table ..."
	
	#data_640K = h5file.getNode('/root').readWhere('(w < 0.792) & (w >= 0.786)')['rg']
	data_602K = h5file.getNode('/root').readWhere('(w < 0.841) & (w > 0.830)')['rg']
	data_300K = h5file.getNode('/root').readWhere('(w > 1.666) &  (w < 1.689)')['rg']
	
	average_300K = numpy.average(data_300K)
	average_602K = numpy.average(data_602K)
	#average_640K = numpy.average(data_640K)

	print data_300K	
	print "<rg> at 300K = ", average_300K
	print data_300K.shape
	print "<rg> at 602K = ", average_602K
	print data_602K.shape
	#print "<rg> at 640K = ", average_640K	
	#print data_640K.shape	

#[Feb/09/2011 10:57:50]        0.836      500000           3     251063.641        602.000           0.808
#[Feb/09/2011 10:57:50]        0.841      500000           3     248190.406        598.000           1.331	
	# [Mar/01/2011 00:37:45]        1.677      500000           3      24297.379        300.000           0.935
	# [Mar/01/2011 00:37:45]        1.689      500000           3      22693.467        298.000           0.923


if __name__ == '__main__':
	main()

