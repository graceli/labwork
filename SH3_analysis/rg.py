import sys
import os
import glob
import copy

import numpy
import tables
import matplotlib

# import my own modules
import plot_and_save2hdf5 as myh5

COUNTER=0
def parse_record(meta, additional):
	"""docstring for parse_record"""
	#chop off record:
	meta_chopped = meta[8:]
	fields = meta_chopped.split("  ")
	# print fields
	# 	print additional
	row = []
	for f in fields:
		key, value = f.split(':')
		row.append(float(value))
	
	rows = []
	newrow = []
	for i in range(0,3):
		newrow = copy.deepcopy(row)
		key, values = additional[i].split(':')
		# print i, values
		data_vector = [ float(v) for v in values.rstrip(' ').lstrip(' ').split(' ') ]
		newrow.extend(data_vector)
		# print i, newrow
		rows.append(newrow)
	
	# for r in rows:
	# 	print r
	# 	
	# sys.exit(0)
	return rows
	
def read_analysis_file(infile):
	"""docstring for read_file"""
	data = []
	while(infile):
		line = infile.readline();
		# print line, len(line)
		if len(line) != 0:
			replica_record = line.rstrip('\n')
			additional_data = []
			for i in range(0,3):
				additional_data.append(infile.readline().rstrip('\n'))
				
			global COUNTER
			print "\rparsing record",COUNTER,
			COUNTER+=1
				
			data.extend(parse_record(replica_record, additional_data))
		else:
			break;
	print
	return data

def parse(datfile, h5file_name):
	"""read all the analysis files into a single h5 file"""
	
	# print "parsing into h5file"
	
	column_names = ['replica', 'sequence', 'w', 'w_nominal', 'rg', 'sas1', 'sas2']
	descr = create_description(column_names, 7)
	h5file = myh5.initialize(h5file_name)
		
	f = open(datfile)
	data = read_analysis_file(f)
	f.close()
	data_array = numpy.array(data)
	myh5.save(h5file, numpy.array(data), '/test', table_struct=descr)

def create_description(column_keys, num_cols, format=tables.Float32Col(dflt=0.0)):
	# print format
	descr = {}
	for i in range(0, num_cols):
		colname = column_keys[i]
		descr[colname] = tables.Float64Col(dflt=0.0, pos=i)
	return descr


def main():
	# Design:
	# parse analyse_force_database files and write into a h5 file
	# for each temperature compute the average and distribution (histogram) for the radius of gyration
	
	temperature_list = [300, 600]
	datfile = sys.argv[1]
	h5_name = sys.argv[2]

	h5file = parse(datfile, h5_name)
	# for T in temperature_list:
	# 	rg = h5file.where(temperature=T)
		# results.append([T, average(rg)])
	# 	print T, average(rg)
	# 	distribution = matplotlib.hist(rg)
	# 	matplotlib.plot(distribution)
	# 
	# matplotlib.savefig('distribution.png')
	
if __name__ == '__main__':
	main()
	
	
