#!/usr/bin/env python
import tables
import os
import sys
import numpy
import pylab

# Example (http://www.pytables.org/docs/manual/ch03.html#id328514)
# Getting object metadata
# >>> print "Object:", table
# Object: /detector/readout (Table(10,)) 'Readout example'
# >>> print "Table name:", table.name
# Table name: readout
# >>> print "Table title:", table.title
# Table title: Readout example
# >>> print "Number of rows in table:", table.nrows
# Number of rows in table: 10
# >>> print "Table variable names with their type and shape:"
# Table variable names with their type and shape:
# >>> for name in table.colnames:
#       print name, ':= %s, %s' % (table.coldtypes[name],
#                                  table.coldtypes[name].shape)
# 
# ADCcount := uint16, ()
# TDCcount := uint8, ()
# energy := float64, ()
# grid_i := int32, ()
# grid_j := int32, ()
# idnumber := int64, ()
# name := |S16, ()
# pressure := float32, ()

def create_description(column_key, num_cols, format=tables.Int32Col(dflt=0)):
	descr = {}
	for i in range(0, num_cols):
		colname = column_key+str(i)
		descr[colname] = format
	
	return descr
	
def save(h5file, data, table_path, table_struct=numpy.dtype(numpy.int64)):
	"""	
		save a numpy array into a given table with name table_name and 
		description table_struct (no compression is used)
	"""
	if type(data) == numpy.ndarray:
		nrows,ncols = data.shape
		table_struct = create_description("col", ncols)
	else:
		print len(data[0])
		print data[0]
		print table_struct.keys()
		print len(table_struct.keys())
		
	# table_path = '/%(group_name)s/%(table_name)s' % vars()
	group_name, table_name = os.path.split(table_path)
	if not h5file.__contains__(table_path):
		if not h5file.__contains__('%(group_name)s' % vars()):
			print "group didn't exist; creating group", group_name
			root, name = os.path.split(group_name)
			group = h5file.createGroup(root, name)
	
		print "table does not exist; creating table", table_path

		table = h5file.createTable(group_name, table_name, table_struct)
	else:
		table = h5file.getNode(table_path)
	
	table.append(data)
	table.flush()
	
	# print "Successfully inserted data with dimensions", data.shape, " into", table_path
			
def initialize(h5_filename, groupName='/'):
	""" 
		open or create a h5 file with predefined
	 	groups inositol, peptide, residue
	"""
	
	h5file = tables.openFile(h5_filename, mode="a")
	filters = tables.Filters(complevel=8, complib='zlib')
	if not h5file.__contains__('/'+groupName):
		h5file.createGroup(h5file.root, groupName, filters=filters)

	return h5file

def main():
	"""
		clean up my flat file analysis and save these flat data files into 
		pytable datafiles
		use ptdump, h5ls tools to read file from commandline
	"""	
	
	if len(sys.argv) < 4:
		print "usage: <file-to-save> <group name> <table name>"
		sys.exit(0)

	file = sys.argv[1]
	groupName = sys.argv[2]
	tableName = sys.argv[3]
	data = numpy.genfromtxt(file, dtype=numpy.float64)

	nrows, ncols = data.shape

	# save our data read in to a *.h5 (HDF5 format)
	# note that name is hard coded in
	analysis_fname = "analysis.h5"
	h5file = initialize(analysis_fname, groupName)
	table_descr = create_description("col", ncols, format=tables.Float64Col(dflt=0.0))

	save(h5file, data, groupName, tableName, table_descr)
	
	# example of how you would read the pytable assuming that you have a table
	# in a file h5file and a group called inositol and a table called inos_bb
	# readout = h5file.root.inositol.inos_bb.read()
	# readout is an numpy array	

if __name__ == "__main__":
	main()

