import glob
import numpy
import sys

isomer = sys.argv[1]
#print isomer
expr = '*%(isomer)s*.dat' % vars()

#print expr
list=glob.glob(expr)
#print list

arraylist = []
for file in list:
	A=numpy.genfromtxt(file,dtype=int)

	# sum over all the elements over axis=0 without the first column
	sum = numpy.sum(A[:,1:],axis=0)

	# change shape so that each amino acid is on a row
	sum.shape = (sum.size/4, 4)
	average_over_peptides = numpy.average(sum,axis=1)
	arraylist.append(average_over_peptides)

# convert to numpy array - needed?
nparray = numpy.array(arraylist)
numpy.savetxt('%(isomer)s_counts.txt' % vars(), nparray, fmt='%0.3f')
print "saved a np array with shape", nparray.shape

# average over all the systems; each system is a row in nparray
average = numpy.average(nparray, axis=0)
print "average", average

std = numpy.std(nparray, axis=0)
print "std", std

#save the normalized average and std
numpy.savetxt('%(isomer)s.txt' % vars(), [(1.0/84700)*average, (1.0/84700)*std], fmt='%0.8f')

