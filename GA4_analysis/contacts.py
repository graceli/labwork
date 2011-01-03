#!/usr/bin/env python

import numpy
import os
import numpy

"""
Contact Analysis for disordered oligomer systems (GA4)
post-processes (averages and standard deviation) data files for polar and nonpolar contacts.
each file is a columnar data file (time series data).
"""
def save(name, A):
	average = numpy.mean(A, axis=1)
	stddev = numpy.std(A, axis=1)

	print name, numpy.mean(A, axis=0), numpy.mean(numpy.mean(A, axis=0)), numpy.std(numpy.mean(A,axis=0))
	numpy.savetxt(name, numpy.transpose(numpy.vstack((average, stddev))))

def nonpolar():
	for system in ['systematic_ap1f', 'systematic_ap2f', 'systematic_4oct']:
		for iso in ['scyllo', 'chiro', 'control']:
			interList = []
			intraList = []
			nrows_min = 200000
			for i in range(0,5):
				inter_filename = '%(system)s_%(iso)s%(i)d_inter_dists_cutoff.dat.1000' % vars()
				intra_filename = '%(system)s_%(iso)s%(i)d_intra_dists_cutoff.dat.1000' % vars()

				if os.path.exists(inter_filename) and os.path.exists(intra_filename):
					print 'analyzing', inter_filename
					print 'analyzing', intra_filename
					inter_data_array = numpy.genfromtxt(inter_filename)
					intra_data_array = numpy.genfromtxt(intra_filename)
					nrows, ncols = inter_data_array.shape
					print nrows, ncols
					if nrows < nrows_min:
						nrows_min = nrows
						
					interList.append(numpy.genfromtxt(inter_filename))
					intraList.append(numpy.genfromtxt(intra_filename))
			
			interList_repack = []
			intraList_repack = []
			if interList:
				for a,b in zip(interList, intraList):
					interList_repack.append(a[0:nrows_min])
					intraList_repack.append(b[0:nrows_min])

			inter_array = numpy.transpose(numpy.array(interList_repack))
			inter_average = numpy.average(inter_array, axis=2)
			print inter_average
			intra_array = numpy.transpose(numpy.array(intraList_repack))
			intra_average = numpy.average(intra_array, axis=2)
			print intra_average
			sum = numpy.transpose(inter_average+intra_average)
			sum[:,0] = sum[:,0]/2
			print sum
			numpy.savetxt('%(system)s_%(iso)s_average.txt' % vars(), sum, fmt='%0.3f')
			

def polar():
	# copied out of systematic_analysis_total.pl
	# $intra[$cols[1]]++;
	# $inter[$cols[2]]++;
	# $total_inos[$cols[3]]++;
	# $per_inositol[$cols[4]]++;
	# $chains_bound[$cols[5]]++;
	# $per_inositol[$cols[6]]++;
	# $chains_bound[$cols[7]]++;
	
	for system in ['systematic_4oct', 'systematic_ap1f', 'systematic_ap2f']:
		for iso in ['scyllo', 'chiro', 'control']:
			data = []
			nrows_min = 200000
			for i in range(0,6):
				filename = '%(system)s_%(iso)s_system%(i)d_all_system_contact_total.moving500' % vars()
				path = os.path.join('data', filename)
				if os.path.exists(path):
					print 'analyzing', filename
					data_array = numpy.genfromtxt(path)
					nrows, ncols = data_array.shape
					if nrows < nrows_min:
						nrows_min = nrows
					data.append(data_array)
			if data:
				# repack list of data array
				repacked = []
				for a in data:
					repacked.append(a[0:nrows_min])

				print numpy.array(repacked).shape
				average = numpy.average(repacked, axis=0)
				stddev = numpy.std(repacked, axis=0)

				numpy.savetxt('%(system)s_%(iso)s_average.txt' % vars(), average, fmt='%0.3f')
				numpy.savetxt('%(system)s_%(iso)s_std.txt' % vars(), stddev,fmt='%0.3f')

if __name__ == '__main__':
	nonpolar()
	
	




