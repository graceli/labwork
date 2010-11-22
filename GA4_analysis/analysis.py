#!/usr/bin/env python

import tables
import numpy

def intersect(polar, nonpolar):
	nrows, ncols = nonpolar.shape
	
	print "nonpolar has dimensions ", nonpolar.shape
	print "polar has dimensions", polar.shape
	print nonpolar
	print polar
	print "polar ", polar[0], polar[0][0]
	
	p_and_np_count  = 0
	nonpolar_count = 0
	polar_count = 0
	total_inositol = 0
	for r in range(0,nrows):
		for i in range(0, ncols):
			# print nonpolar_new[r][i], polar[r][i]
			if nonpolar[r][i] > 0.0 and polar[r][i] > 0.0:
				p_and_np_count += 1
			elif nonpolar[r][i] > 0.0 and polar[r][i] == 0.0:
				nonpolar_count += 1
			elif nonpolar[r][i] == 0.0 and polar[r][i] > 0.0:
				polar_count += 1
			total_inositol += 1
	
	return polar_count, nonpolar_count, p_and_np_count, total_inositol

def binding(polar, nonpolar):
	nrows, ncols = nonpolar.shape
	
	total_binding = nonpolar + polar
	
	# aggregate is a single column of 0s or non-zeros indicating bound and unbound
	aggregate = numpy.sum(total_binding, axis=1)
	nonzero_indices = numpy.nonzero(aggregate)
	
	# bound is the number of nonzero numbers in the array of bound and unbound
	bound = len(aggregate[nonzero_indices])
	unbound = len(aggregate) - bound
	
	return bound, unbound
	
			
def getTable(h5file,path):
	if h5file.__contains__(path):
		return h5file.getNode(path)
	else:
		print path, " table does not exist in h5file"
		return None

def convert_to_numpy(table):
	numpy_array = table.read().view(dtype=numpy.float64).reshape(-1, len(table[0]))
	return numpy_array
	

def disordered():
	"""This is computes the intersection of polar and nonpolar binding of inositol"""
	
	polar_h5 = tables.openFile('GA4_disordered_polar_analysis.h5', mode='a')
	nonpolar_h5 = tables.openFile('GA4_disordered_nonpolar_analysis.h5', mode='a')

	f = open('data.txt', 'w')
	print >>f, "#table_name polar nonpolar p_and_np polar_fraction nonpolar_fraction p_and_np_fraction total_inositols bound unbound"
	for system in ["ap1f", "oct"]:
		for iso in ["scyllo", "chiro"]:
			for i in range(0,5):
				table_path = '/%(system)s/%(iso)s%(i)d' % vars()
				polar_table = getTable(polar_h5, table_path)
				nonpolar_table = getTable(nonpolar_h5, table_path)
				
				if polar_table != None and nonpolar_table != None:
					polar_array = convert_to_numpy(polar_table)
					nonpolar_array = convert_to_numpy(nonpolar_table)
					polar,nonpolar,p_and_np, total_inositol = intersect(polar_array, nonpolar_array)
					bound,unbound = binding(polar_array, nonpolar_array)
					total = float(polar + nonpolar + p_and_np)
					print >>f, table_path, polar, nonpolar, p_and_np, polar/total, nonpolar/total, p_and_np/total, total_inositol, bound, unbound, unbound/float(bound)*125 

def mon():
	"""This is computes the intersection of polar and nonpolar binding of inositol"""
	
	polar_h5 = tables.openFile('GA4_mon_polar_analysis.h5', mode='a')
	nonpolar_h5 = tables.openFile('GA4_mon_nonpolar_analysis.h5', mode='a')
	# f = open('data.txt', 'w')
	# print >>f, "#table_name polar nonpolar p_and_np polar_fraction nonpolar_fraction p_and_np_fraction total_inositols bound unbound K"
	for system in ["mon"]:
		for iso in ["scyllo", "chiro"]:
			data = []
			for i in range(0,1118):
				table_path = '/%(system)s/%(iso)s%(i)d' % vars()
				print "analyzing",table_path
				polar_table = getTable(polar_h5, table_path)
				nonpolar_table = getTable(nonpolar_h5, table_path)
				if polar_table != None and nonpolar_table != None:
					polar_array = convert_to_numpy(polar_table)
					nonpolar_array = convert_to_numpy(nonpolar_table)
					print "computing intersection ..."
					polar, nonpolar, p_and_np, total_inositol = intersect(polar_array[0:5001,1:], nonpolar_array[0:5001,1:])
					
					print "computing binding constant ..."
					bound,unbound = binding(polar_array[0:5001,1:], nonpolar_array[0:5001,1:])
					total = float(polar + nonpolar + p_and_np)
					# print >>f, table_path, polar, nonpolar, p_and_np, polar/total, nonpolar/total, p_and_np/total, total_inositol, bound, unbound, unbound/float(bound)*125
					data.append(numpy.array([polar, nonpolar, p_and_np, polar/total, nonpolar/total, p_and_np/total, total_inositol, bound, unbound, unbound/float(bound)]))
				else:
					print "table not found"
		
			data = numpy.vstack(data)
			numpy.savetxt('%(iso)s_data.txt' % vars(), data, fmt='%-0.3f')
			data_aggregation = numpy.sum(data, axis=0)
			numpy.savetxt('%(iso)s_data_aggregation.txt' % vars(), numpy.transpose(data_aggregation), fmt='%.3f')
	
def main():
	mon()
	
if __name__ == '__main__':
	main()
