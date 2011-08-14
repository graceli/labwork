import tables
import os
import sys
import numpy

def create_description(column_key, num_cols, format=tables.Int32Col(dflt=0)):
	descr = {}
	# descr.update(metacols)
	for i in range(0, num_cols):
		colname = column_key+str(i)
		descr[colname] = tables.Int32Col(dflt=0.0, pos=i)
	
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
	filters = tables.Filters(complevel=8, complib='zlib')
	for group_name in ["inositol", "peptide", "residue"]:
		if not h5file.__contains__('/' + group_name):
			h5file.createGroup(h5file.root, group_name, filters=filters)

	return h5file

def main():
	"""
		clean up my inositol flat file analysis and save these flat data files into 
		pytable and HDF5 files

	"""	
	if len(sys.argv) < 2:
		print "Error: missing input files"
		sys.exit(1)
		
	path = sys.argv[2]
	analysis_file = sys.argv[1]
	isomer = sys.argv[3]

	# tables_names = ["inositol-inos_bb", "inositol-inos_glu", "inositol-inos_lys", "inositol-""peptide-p2p_vs_t", "peptide-pep_bb", "peptide-pep_side", "residue-res_bb", "residue-res_side", "residue-per_res_contacts", "residue-per_inos_contacts"]
	tables_names = ["inositol-inos_total", "residue-per_inos_contacts"]
	#group_name = {'inos': 'inositol', 'pep' : 'peptide', 'res' : 'residue'}
	table_descr = {}
	
	h5file = initialize(analysis_file)
	
	for name in tables_names:
		for i in range(1,600):
			#file_base = "sys%(i)s_nosol.xtc_" % vars()
			file_base = isomer + "_sys%(i)s_" % vars()
			#print name.split("-")
			group_name,extension = name.split("-")
			filename = file_base + extension + ".dat"
			col_key = extension.split("_")[0]
			
			path_to_file = os.path.join(path,filename)
			if os.path.exists(path_to_file):
				print "loading", path_to_file
				data = numpy.genfromtxt(path_to_file, dtype=int, comments="#")[0:7500]
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
				save(h5file, data, group_name, extension, table_descr[extension])

	# readout = h5file.root.inositol.inos_bb.read()
	# print numpy.asscalar(numpy.array(readout[0][0]))
	
if __name__ == '__main__':
	main()
