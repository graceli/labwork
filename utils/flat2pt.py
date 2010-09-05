import tables
import os
import sys
import numpy

def create_description(column_key, num_cols, format=tables.Int32Col(dflt=0)):
	descr = {}
	# descr.update(metacols)
	print descr
	for i in range(0, num_cols):
		colname = column_key+str(i)
		descr[colname] = format
	
	return descr
	
def save(h5file, data, group_name, table_name, table_struct):
	"""	
		save a numpy array into a given table with name table_name and 
		description table_struct
	"""

	table_path = '/%(group_name)s/%(table_name)s' % vars()
	
	# try:
	# 	table = h5file.getNode(table_path)
	# except NoSuchNodeError:
	# 	print "table does not exist; creating table", table_name
	# 	table = h5file.createTable('/' + group_name, table_name, table_struct)
			
	if not h5file.__contains__(table_path):
		print "table does not exist; creating table", table_name
		table = h5file.createTable('/' + group_name, table_name, table_struct)
	else:
		table = h5file.getNode(table_path)
	
	table.append(data)
	table.flush()
	
	print "inserted data with dimensions", data.shape, " into", table_path
			
def initialize(h5_filename):
	""" 
		open or create a h5 file with predefined
	 	groups inositol, peptide, residue
	"""
	
	h5file = tables.openFile(h5_filename, mode="a")
	for group_name in ["inositol", "peptide", "residue"]:
		if not h5file.__contains__('/' + group_name):
			h5file.createGroup(h5file.root, group_name)

	return h5file

def main():
	"""
		clean up my inositol flat file analysis and save these flat data files into 
		pytable and HDF5 files

	"""	
	
	path = sys.argv[1]
	
	print "path", path
	
	tables_names = ["inos_bb", "inos_glu", "inos_lys", "pep_p2p_vs_t", "pep_bb", "pep_side", "res_bb", "res_side"]
	group_name = {'inos': 'inositol', 'pep' : 'peptide', 'res' : 'residue'}
	table_descr = {}
	
	h5file = initialize('analysis.h5')
	
	for extension in tables_names:
		for i in range(1,500):
			file_base = "sys%(i)s_nosol.xtc_" % vars()
			# print "load", file_base + extension + ".dat"
			filename = file_base + extension + ".dat"
			col_key = extension.split("_")[0]
			
			path_to_file = os.path.join(path,filename)
			if os.path.exists(path_to_file):
				print "loading", path_to_file
				data = numpy.genfromtxt(path_to_file, dtype=int, comments="#")
				print "loaded a numpy array with shape", data.shape
				shape = data.shape
				
				if data.ndim == 1:
					nrows = data.size
					col_num = 1
				else:
					nrows = shape[0]
					col_num = shape[1]
				
				# table_descr[extension] = create_description(col_key, col_num+1 ,{'sys':tables.Int32Col(dflt=0)})
				table_descr[extension] = create_description(col_key, col_num)
				save(h5file, data, group_name[col_key], extension, table_descr[extension])

	# readout = h5file.root.inositol.inos_bb.read()
	# print numpy.asscalar(numpy.array(readout[0][0]))
	
if __name__ == '__main__':
	main()
	
	