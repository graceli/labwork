#!/usr/bin/env python
import tables
import os
import sys
import numpy
import ConfigParser
import logging

def create_description(column_key, num_cols, format='int'):
	descr = {}
	# descr.update(metacols)
	for i in range(0, num_cols):
		colname = column_key+str(i)
		descr[colname] = get_data_type(format, col_pos=i)	
	return descr
	
def save(h5file, data, group_name, table_name, table_struct=None):
	"""	
	save a numpy array into a given table with name table_name and 
	description table_struct
	"""

	# try:
	# 	table = h5file.getNode(table_path)
	# except NoSuchNodeError:
	# 	print "table does not exist; creating table", table_name
	# 	table = h5file.createTable('/' + group_name, table_name, table_struct)

	group_path = '/%(group_name)s' % vars()
	if not h5file.__contains__(group_path):
		logging.info("group does not exist; creating group: %s", group_name)
		h5file.createGroup(h5file.root, group_name)

	table_path = '/%(group_name)s/%(table_name)s' % vars()
	if not h5file.__contains__(table_path):
		logging.info("table does not exist; creating table %s", table_name)
		#setting the description with an numpy array creates a description with the numpy array dtype
		#and injects the data in the array into the table, if any
		# print data.dtype
		logging.info("creating Table")
		table = h5file.createTable('/' + group_name, table_name, description=table_struct)
		# table = h5file.createArray('/'+group_name, table_name, data)
	else:
		table = h5file.getNode(table_path)
	
	table.append(data)
	table.flush()
	
	logging.info("inserted data with dimensions %s into %s", data.shape, table_path)
			
def initialize(h5_filename, groups=[]):
	""" 
		open or create a h5 file with predefined
	 	groups inositol, peptide, residue
	"""
	filters = tables.Filters(complevel=8, complib='zlib')
	h5file = tables.openFile(h5_filename, mode="a", title='analysis', filters=filters)

	return h5file

def get_data_type(type, col_pos=0):
	data_types = {'int': tables.Int64Col(dflt=0, pos=col_pos), 'float': tables.Float64Col(dflt=0.0, pos=col_pos)}
	return data_types[type]

def numpy_dtype(type):
	data_types = {'int': numpy.int64, 'float': numpy.float64 }
	return data_types[type]
	
def main():
	
	filename = sys.argv[1]
	
	LOG_FILENAME = filename+'.log'

	
	logging.basicConfig(level=logging.INFO,
	                    format='%(asctime)s %(levelname)-8s %(message)s',
	                    datefmt='%a, %d %b %Y %H:%M:%S',
	                    filename=LOG_FILENAME,
	                    filemode='w')
	
	logging.info(" ==================== starting abeta_analysis.py ====================")
	logging.info(" Data will be written to %s with messages in %s", filename, LOG_FILENAME)
	
	config = ConfigParser.ConfigParser()
	base,name=os.path.split(filename)
	read = config.read(os.path.join(base,'config.ini'))
	
	h5file = initialize(filename)
	sections = config.sections()

	groups = sections
	groups.remove('analysis')
	
	logging.info("Found groups %s", groups)
	for analysis_section, value in config.items('analysis'):
		# print analysis_section, value
		if value == "False":
			logging.info("%s set to false. Removing %s", analysis_section, analysis_section)
			groups.remove(analysis_section)
		
	logging.info("Groups %s will be analyzed", groups)
	for group_name in groups:
		logging.info("reading files under section %s", group_name)
		fmt=config.get(group_name, 'dtype')
		config.remove_option(group_name, 'dtype')
		# print dtype, data_types[dtype]
		# print config.items(group_name)
		# sys.exit(0)
		for table_name, files in config.items(group_name):
			l = files.split(' ')
			print l
			
			all_data=[]
			all_data_matrix=None
			rows=0
			cols=0
			if l == None or len(l) == 1:
				logging.info("saving %s in %s %s", files, group_name, table_name)
				all_data_matrix = numpy.genfromtxt(files)
				# rows,cols = all_data_matrix.shape
				# 		logging.info("shape: %s", shape)
			else:
				logging.info("found %d files to read and save", len(l))
				for datafile in l:
					logging.info("saving %s in %s %s", datafile, group_name, table_name)
					data = numpy.genfromtxt(datafile)
					all_data.append(data)
					all_data_matrix = numpy.hstack(all_data)
			
			rows,cols = all_data_matrix.shape
			descr = create_description('col', cols, format=fmt)
			logging.info("description created %s", descr)
			save(h5file, all_data_matrix, group_name, table_name, descr)
	
	logging.info("finished processing files")

if __name__ == '__main__':
	main()


	
	# ro = file.getNode('/ab_15_scyllo/inositol_residue_np_ch5').read()
	# ro.view(dtype=np.int64).reshape(-1,len(ro[0]))
	# sum(ro.view(dtype=np.int64).reshape(-1,len(ro[0])),axis=0)
	
	#organization of the h5 file
	#scyllo (group)
	#table
	#protein
	#col=rmsd col=SASA
	#inositol_polar
	#col = inositol contact numbers
	#inositol_nonpolar
	#col = inositol contact numbers
	#chiro (group)
	#water (group)
	#protein
	
	# algorithm:
	# given a single XTC and a set of system parameters
	# do gromacs analysis on the dataset
	# parse and process gromacs analysis text files
	# tar and gz the analysis text files
	# write data to h5 file
	#	use compression
	# use proper logging and config parsing
	# idea write a base class from which I extend each time I create a new analysis
	# how to override base class methods
	# how to create packages




