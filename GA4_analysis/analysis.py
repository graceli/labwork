#!/usr/bin/env python

import tables
import numpy
import math
import plot_and_save2hdf5 as myh5
import re
import csv
import glob
import os

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
	print "number of rows", nrows
	# total_binding = nonpolar + polar
	
	# aggregate is a single column of 0s or non-zeros indicating bound and unbound
	aggregate = numpy.sum(polar+nonpolar, axis=1)
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

def convert_to_numpy(table, dtype=numpy.float64):
	numpy_array = table.read().view(dtype=dtype).reshape(-1, len(table[0]))
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

					
def binding_error(polar_array, nonpolar_array, block_size=1000):
	nrows,ncols = polar_array.shape
	max_num_blocks = int(math.floor(nrows/block_size))
	
	Kd_blocklist = []
	for i in range(0,max_num_blocks):
		block_start = block_size*i
		block_end = min(5001, block_size*(i+1)+1)
		print "block_start=", block_start, "block_end=", block_end
		print polar_array[block_start:block_end, 1:]
		print numpy.array(binding(polar_array[block_start:block_end, 1:], nonpolar_array[block_start:block_end, 1:]))
		bound_unbound_vector = numpy.array(binding(polar_array[block_start:block_end, 1:], nonpolar_array[block_start:block_end, 1:]))
		Kd_blocklist.append(bound_unbound_vector)

	return numpy.hstack(Kd_blocklist)

def block_average(data):
	print "computing block averaging"
	print "the data has shape", data.shape
	block_size = 200
	nblocks = int(math.floor(1117/block_size))
	blocklist = []
	for i in range(0, nblocks):
		start = i*block_size
		end = (i+1)*block_size
		block_average = numpy.sum(data[start:end, :], axis=0)
		blocklist.append(block_average)
	
	return numpy.vstack(blocklist)
							
def mon(offset=0, max_num_dataset=10):
	"""This is computes the intersection of polar and nonpolar binding of inositol"""
	
	polar_h5 = tables.openFile('GA4_mon_polar_analysis.h5', mode='a')
	nonpolar_h5 = tables.openFile('GA4_mon_nonpolar_analysis.h5', mode='a')

	# print >>f, "#table_name polar nonpolar p_and_np polar_fraction nonpolar_fraction p_and_np_fraction total_inositols bound unbound K"

	for system in ["mon"]:
		for iso in ["scyllo", "chiro"]:
			data = []
			errorlist = []
			for i in range(0, max_num_dataset):
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
					bound,unbound = binding(polar_array[offset:5001, 1:], nonpolar_array[offset:5001, 1:])	
					# errorlist.append(binding_error(polar_array[offset:5001, 1:], nonpolar_array[offset:5001, 1:], block_size=1000))
					
					total = float(polar + nonpolar + p_and_np)
					# print >>f, table_path, polar, nonpolar, p_and_np, polar/total, nonpolar/total, p_and_np/total, total_inositol, bound, unbound, unbound/float(bound)*125
					data.append(numpy.array([polar, nonpolar, p_and_np, total, total_inositol, bound, unbound]))
				else:
					print "table not found"
		
			data = numpy.vstack(data)
			# numpy.savetxt('%(iso)s_data.txt' % vars(), data, fmt='%-0.3f')
			data_sum = numpy.sum(data, axis=0)
			
			# nasty output
			heading = ["polar", "nonpolar", "p_and_np", "total", "total_inositol", "nbound", "nunbound"]
			d = dict(zip(heading, data_sum))
			f = open('%(iso)s_data_aggregation.txt' % vars(), 'w')
			for key in sorted(d.keys()):
				print >>f, key, d[key]
			f.close()
			
			#compute error
			blocklist = block_average(data)
			numpy.savetxt('%(iso)s_blocks.txt' % vars(), blocklist, fmt='%.3f')
						
			# numpy.savetxt('%(iso)s_data_aggregation.txt' % vars(), numpy.transpose(data_aggregation), fmt='%.3f')
			#output the error
			# error_array = numpy.vstack(errorlist)
			# error_aggregate = numpy.sum(error_array, axis=0)
			# numpy.savetxt('%(iso)s_data_aggregation_Kd_error.txt' % vars(), numpy.vstack(error_aggregate), fmt='%d')	

def stoichiometry(polar_array, nonpolar_array):
	""" computes the binding stoichiometry of inositol
		ie. reports the number of molecules bound to the aggregate/protein
	"""
	print polar_array.shape
	print nonpolar_array.shape

	total = polar_array + nonpolar_array

	nrows, ncols = total.shape
	counts = [0,0,0]
	for i in range(0, nrows):
		if (total[i][0] and total[i][1]):
			counts[2] += 1
		elif (not total[i][0] and not total[i][1]):
			counts[0] += 1
		else:
			counts[1] += 1
	return counts

def analysis(saveto_h5, max_num_dataset=10):
	""" A bad way to organize a sequence of analysis """
	
	# h5 files to read, tables, and paths to tables are encoded inside the analysis
	# ideally they would be refactored into a configuration file
	polar_h5 = tables.openFile('GA4_mon_polar_analysis.h5', mode='a')
	nonpolar_h5 = tables.openFile('GA4_mon_nonpolar_analysis.h5', mode='a')

	# analyze and aggregate all data for each iso and store each in a separate table 
	for system in ["mon"]:
		for iso in ["scyllo", "chiro"]:
			# clear the results
			analysis_results = []
			for i in range(0, max_num_dataset):
				table_path = '/%(system)s/%(iso)s%(i)d' % vars()
				print "analyzing", table_path
				polar_table = getTable(polar_h5, table_path)
				nonpolar_table = getTable(nonpolar_h5, table_path)
				if polar_table != None and nonpolar_table != None:
					polar_array = convert_to_numpy(polar_table)
					nonpolar_array = convert_to_numpy(nonpolar_table)
					s = stoichiometry(polar_array[0:5001, 1:], nonpolar_array[0:5001,1:])
					analysis_results.append(s)

			myh5.save(saveto_h5, numpy.vstack(analysis_results), "/mon_analysis/stoichiometry_%(iso)s" % vars())

def aggregate(h5, aggregate_function=numpy.sum, where='/'):
	# sum over all rows of each table in the file
	# and print to screen the result
	for t in h5.listNodes(where):
		a = numpy.sum(convert_to_numpy(t, dtype=numpy.int32), axis=0)
		sum_across = aggregate_function(a)
		print a/float(sum_across)

def distribution(h5, iso="scyllo", where='/'):
	total_counts = [0] * 9
	total_points = 0
	num_tables = 0

	pattern = re.compile(r'%(iso)s' % vars())
	for table in h5.listNodes(where):
		# only process the table that matches the isomer
		# this is kind of hacky
		match = pattern.search(table._v_name)
		if match:
			array = convert_to_numpy(table)
			no_time = array[:,1:]
			nrows, ncols = no_time.shape
			for col in range(0, ncols):
				column = no_time[:, col]
				counts, bins = numpy.histogram(column, bins=range(0,10))
				total_counts += counts
				total_points += 5001		
			num_tables += 1
	print "total number of tables processed", num_tables
	numpy.savetxt('%(iso)s_distribution.txt' % vars(), numpy.transpose(numpy.vstack([range(0,9), total_counts/float(total_points)])), fmt='%0.3f')
	
def dssp():
	# "Coil": Float64Col(shape=(), dflt=0.0, pos=2),
	#   "B-Sheet": Float64Col(shape=(), dflt=0.0, pos=3),
	#   "B-Bridge": Float64Col(shape=(), dflt=0.0, pos=4),
	#   "Bend": Float64Col(shape=(), dflt=0.0, pos=5),
	#   "Turn": Float64Col(shape=(), dflt=0.0, pos=6),
	#   "A-Helix": Float64Col(shape=(), dflt=0.0, pos=7),
	#   "3-Helix": Float64Col(shape=(), dflt=0.0, pos=8),
	
	""" 
		reads the h5 file containing the secondary structure analysis timeseries
		and outputs a distribution for each datafile
		write the computed distribution to a csv file with headings
	"""
	f = tables.openFile('analysis_results.h5')
	structList = ['Coil', 'B-Sheet', 'B-Bridge', 'Bend', 'Turn', 'A-Helix', '3-Helix']
	for system in ['4oct', 'ap1f', 'ap2f']:
		for iso in ['scyllo', 'chiro', 'control']:
			pattern = re.compile(r"%(system)s_%(iso)s" % vars())
			ss_distribution = {'#filename': '', 'Coil' : 0.0, 'B-Sheet' : 0.0, 'B-Bridge' : 0.0, 'Bend' : 0.0, 'Turn' : 0.0, 'A-Helix' : 0.0, '3-Helix' : 0.0}

			# use regex to match the tables for 'system' and 'iso' -- another better way? 
			heading = ['#filename', 'Coil', 'B-Sheet', 'B-Bridge', 'Bend', 'Turn', 'A-Helix', '3-Helix']
		
			#write data to ascii csv to import into Excel
			dictFilename = '%(system)s_%(iso)s_dict.csv' % vars()
			fhandle = open(dictFilename, 'wb')
		 	writer = csv.DictWriter(fhandle, heading, delimiter=" ")
			writer.writerow(dict(zip(heading, heading)))

			num_matches = 0
			for t in f.listNodes(where='/dssp'):
				filename = t.col('filename')[0]
				match = pattern.search(filename)
				if match:
					ss_distribution['#filename'] = filename
					print "processing", filename
					# compute the sum of ss distribution over
					# all independent runs
					for struct in ss_distribution.keys():
						if struct != '#filename':
							try:
								ss = t.col(struct)[0]
							except KeyError:
								print filename, struct,"does not exist in table"
							else:
								# print filename, struct, ss
								ss_distribution[struct] = float(ss)
						else:
							continue
					writer.writerow(ss_distribution)
					fhandle.flush()
			fhandle.close()
			#regenerate csv file and compute the average and std (over all independent systems)
			#TODO: figure out how to do this without read off of disk
			#get rid of the first column because its all strings (read in as NaN by Numpy)
	
	data = []
	for dictFilename in glob.glob("*.csv"):
		try:
			data = numpy.genfromtxt(dictFilename, comments='#')[:,1:]
		except IOError:
			print data
			print dictFilename, "size", os.path.getsize(dictFilename), "had a problem"


		# print data
		average = numpy.average(data, axis=0)
		print average
		stdev = numpy.std(data, axis=0)
		writer = csv.DictWriter(open(dictFilename, 'a'), heading, delimiter=" ")
		#append the average and stdev back into the csv file
		aggregate_data = {'average':average, 'stdev':stdev}
		for newhead in ['average', 'stdev']:
			index = 0
			for attr in heading:
				if attr == "#filename":
					ss_distribution[attr] = newhead
				else:
					print index
					ss_distribution[attr] = aggregate_data[newhead][index]
					index += 1
			writer.writerow(ss_distribution)

def main():
	# h5file = tables.openFile('GA4_mon_polar_analysis.h5')
	# mon(max_num_dataset=1117)
	# saveto_h5 = myh5.initialize('analysis_results.h5')
	# analysis(saveto_h5, max_num_dataset=1117)
	# aggregate(h5file, where='/mon_analysis')
	# distribution(h5file, iso="chiro", where='/mon/')
	# distribution(h5file, iso="scyllo", where='/mon/')
	
	dssp()

if __name__ == '__main__':
	main()
