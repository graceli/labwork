#!/usr/bin/env python
import pylab
import numpy
import csv
import re
import glob
import tables
def build_plot_data(h5file, tables):
	""" Builds a list of data for each given table to be plotted
		This changes depending on the analysis """
	
	data_all = []
	for table in tables:
		datatable = h5file.getNode(table).read()
		print datatable.dtype
		row_size = len(datatable[0])
		print "build_plot_data:", row_size
		data_matrix = datatable.view(dtype=numpy.int64).reshape(-1, row_size)
		print data_matrix
		nrows, ncols = data_matrix.shape
		# nonpolar_contact_sum.append(list(data.sum(axis=0)[1:]/float(nrows)))
		data_sum = numpy.sum(data_matrix, axis=0)[1:]/float(nrows)  #sum over rows of the matrix, and get rid of the first column 
		data_all.append(data_sum)
	
	return data_all

def plot_matshow_axis():
	# nonpolar_dict={}
	# 		# build dictionary 
	# 		for i in range(0, len(residue_map)):
	# 			print residue_map[i], nonpolar_contact_sum[i]
	# 			nonpolar_dict[residue_map[i]] = nonpolar_contact_sum[i]
	
	
	# #print dictionary in sorted order
	# 		nonpolar_contact_sum_in_order = []
	# 		residue_map_in_order = []
	# 		for key in sorted(nonpolar_dict.keys(),key=lambda k: int(re.search('[0-9]+', k).group(0))):
	# 			print key, nonpolar_dict[key]
	# 			residue_map_in_order.append(key)
	# 			nonpolar_contact_sum_in_order.append(float(nonpolar_dict[key])/int(residue_map_dict[key]))
	pass
	
def plot_matshow(filename, data, data_axis=None):
	""" plots a matrix using the matshow function with a specific formating """

	pylab.rcParams['xtick.labelsize']='8'
	pylab.rcParams['legend.fontsize']='8'	
	pylab.axis('tight')
	
	# pylab.matshow(numpy.transpose(data), cmap=pylab.cm.spectral)
	pylab.matshow(data, cmap=pylab.cm.spectral)
	pylab.colorbar(orientation="horizontal", shrink=0.75, extend='max')
	#pylab.yticks(x, tuple(data_axis))
	pylab.savefig(filename + '.pdf' % vars())	
	numpy.savetxt(filename + '.dat', numpy.transpose(data), fmt='%f')

def main():
	"""docstring for main"""
	
	# systems = ['ab_15_scyllo']
	h5list = glob.glob("*.h5")
	tables_list = ['/ab_15_scyllo/inositol_residue_np_ch1','/ab_15_scyllo/inositol_residue_np_ch2',
	'/ab_15_scyllo/inositol_residue_np_ch3','/ab_15_scyllo/inositol_residue_np_ch4',
	'/ab_15_scyllo/inositol_residue_np_ch5']
	
	print "list of h5 files to read:", h5list
	print "reading tables:", tables_list
	
	all_systems_nonpolar_contact = []
	for h5file in h5list:
		f = tables.openFile(h5file, mode='a')
		data = build_plot_data(f, tables_list)
		all_systems_nonpolar_contact.append(data)

	
	# nonpolar_axis = plot_matshow_axis()
	# 	plot_matshow('test', all_systems_nonpolar_contact, nonpolar_axis)
	
if __name__ == '__main__':
	# systems = ['ab_15_scyllo']
	h5list = glob.glob("*.h5")
	tables_list = ['/ab_15_scyllo/inositol_residue_np_ch1','/ab_15_scyllo/inositol_residue_np_ch2',
	'/ab_15_scyllo/inositol_residue_np_ch3','/ab_15_scyllo/inositol_residue_np_ch4',
	'/ab_15_scyllo/inositol_residue_np_ch5']
	
	print "list of h5 files to read:", h5list
	print "reading tables:", tables_list
	
	all_systems_nonpolar_contact = []
	for h5file in h5list:
		f = tables.openFile(h5file, mode='a')
		data = build_plot_data(f, tables_list)
		all_systems_nonpolar_contact.append(data)

	
	all_data_as_matrix = numpy.array(all_systems_nonpolar_contact)
	scyllo_chain_np_contact = numpy.average(all_data_as_matrix[0:10], axis=0)
	chiro_chain_np_contact = numpy.average(all_data_as_matrix[11:], axis=0)
	#nonpolar_axis = plot_matshow_axis()
	plot_matshow('test_sc', scyllo_chain_np_contact)
	plot_matshow('test_ch', chiro_chain_np_contact)
