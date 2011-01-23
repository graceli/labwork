import plot_and_save2hdf5 as myh5
import glob
import os
import numpy
import tables
import utils

""" Analysis for beta-aggregate for GA4 """

# TODO: Geometric binding mode analysis
# TODO: Frequency of binding to C=O vs N-H (I have to parse this data)
def intersect(polar, nonpolar):
	nrows, ncols = polar.shape
	
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
		# first column is ignored; it is assumed to be time
		for i in range(1, ncols):
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
	print "number of rows", nrows
	# total_binding = nonpolar + polar
	
	# aggregate is a single column of 0s or non-zeros indicating bound and unbound
	# ignore the first column; it is assumed to be time
	try:
		aggregate = numpy.sum(polar[:,1:] + nonpolar[:, 1:], axis=1)
	except ValueError:
		print "polar", polar.shape
		print "nonpolar", nonpolar.shape
	else:
		nonzero_indices = numpy.nonzero(aggregate)
	
	# bound is the number of nonzero numbers in the array of bound and unbound
	bound = len(aggregate[nonzero_indices])
	unbound = len(aggregate) - bound
	
	return bound, unbound

def stoichiometry(polar_array, nonpolar_array):
	""" computes the binding stoichiometry of inositol
		ie. reports the number of molecules bound to the aggregate/protein
	"""
	print polar_array.shape
	print nonpolar_array.shape

	total = polar_array + nonpolar_array
	nrows, ncols = total.shape
	counts = [0]*ncols
	for i in range(0, nrows):
		nbound = 0
		# ignore the first column as time
		for n in range(1, ncols):
			if total[i][n]:
				nbound += 1
		try:
			counts[nbound] += 1
		except IndexError:
			print nbound

	try:
		normalize = counts/numpy.sum(counts)
	except:
		print normalize
		return counts
	else:
		return counts

def read_polar(h5file):
	""" reads in the flat files contain polar contact analysis into a
		h5 file
	"""
	# Contact Analysis (per inositol)
	# Polar contact
	polar = glob.glob("polar/*.dat")
	for file in polar:
		# print file
		data = numpy.genfromtxt(file, dtype=numpy.int32)
		# construct table path from filename
		discard, rest = os.path.split(file)
		parts = rest.split('_')
		group_name, ext = parts[-1].split('.')
		table_name = '_'.join(parts[0:4])
		table_path = os.path.join('/', os.path.join(group_name, table_name))
		
		print "saving %(file)s to" % vars(), table_path
		
		myh5.save(h5file, data, table_path)


def read_nonpolar(h5file):
	""" reads in the flat files containing nonpolar contact analysis into a 
		h5 file
	"""
	nonpolar = glob.glob("nonpolar_all/*per_inositol_contacts.dat")
	for file in nonpolar:
		# print file
		data = numpy.genfromtxt(file, dtype=numpy.int32)
		# Nasty file name string parsing
		path, filename = os.path.split(file)
		parts = filename.split('_')
		sys_number = parts[2][-1]
		table_name = 'inf_' + '_'.join(parts[1:3])[:-1] + '_sys' + sys_number
		group_name = 'nonpolar_per_inositol'
		table_path = os.path.join(os.path.join('/', group_name), table_name)
		print "saving %(file)s to" % vars(), table_path
		myh5.save(h5file, data, table_path)
		
def load():
	# initialize a h5 file to store all the analysis relating to GA4-beta protofibrils
	h5file = myh5.initialize('GA4_beta_analysis.h5')
	read_polar(h5file)
	read_nonpolar(h5file)
	h5file.close()

def analyze(h5file):
	"""This is computes the intersection of polar and nonpolar binding of inositol"""
	# open the h5 file for protofibrils
	inositol_concentration_mM = 119
	f = open('data.txt', 'w')
	print >>f, "#table_name polar nonpolar p_and_np polar_fraction nonpolar_fraction p_and_np_fraction total_inositols bound unbound Kd stoichiometry"
	 
	for t in h5file.listNodes('/nonpolar_per_inositol'):
		# find matching polar table
		nonpolar_table = t
		polar_table = h5file.getNode(os.path.join('/polar-contacts-per-inositol', t.name))
		
		if polar_table != None and nonpolar_table != None:
			polar_array = utils.convert_to_numpy(polar_table)
			nonpolar_array = utils.convert_to_numpy(nonpolar_table)
			polar,nonpolar,p_and_np, total_inositol = intersect(polar_array, nonpolar_array)
			bound,unbound = binding(polar_array, nonpolar_array)
			stoic = stoichiometry(polar_array, nonpolar_array)
			total = float(polar + nonpolar + p_and_np)
			print >>f, polar_table.name, nonpolar_table.name, polar, nonpolar, p_and_np, polar/total, nonpolar/total, p_and_np/total, total_inositol, bound, unbound, unbound/float(bound)*inositol_concentration_mM, stoic

def main():
	load()
	h5file = tables.openFile('GA4_beta_analysis.h5')
	analyze(h5file)
	
if __name__ == '__main__':
	main()

