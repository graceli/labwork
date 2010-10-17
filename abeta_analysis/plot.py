#!/usr/bin/env python
import pylab
import numpy
import csv
import re
import glob
import tables
import sys

def filter_nonpolar(data_matrix):
	return data_matrix

def filter_hbond(data_matrix):
	return data_matrix[:,1::3]
	
def reorder(data, chain_num):
	""" reorder the data according to a textfile mapping """
	#load in residues
	r=csv.reader(open('chain%(chain_num)d_table.dat' % vars()), delimiter=' ')
	residue_map=[]
	residue_map_dict = {}
	for line in r:
		residue_map.append(line[0])
		residue_map_dict[line[0]]=int(line[1])

	nonpolar_dict={}
	# build dictionary 
	for i in range(0, len(residue_map)):
		nonpolar_dict[residue_map[i]] = data[i]

	data_in_order = []
	for key in sorted(nonpolar_dict.keys(),key=lambda k: int(re.search('[0-9]+', k).group(0))):
		#return data in the order of residue number and normalize by the number of atoms 
		data_in_order.append(float(nonpolar_dict[key])/int(residue_map_dict[key]))

	return data_in_order

def build_plot_data(h5file, tables,filter_func=filter_nonpolar):
	""" Builds a list of data for each given table to be plotted
		This changes depending on the analysis """
	
	#read chain tables and reorder data according to the residue mapping	
	data_all = []
	chain_num=1
	for table in tables:
		print table
		datatable = h5file.getNode(table).read()
		# print datatable.dtype
		row_size = len(datatable[0])
		# print "build_plot_data:", row_size
		data_matrix = datatable.view(dtype=numpy.int64).reshape(-1, row_size)
		# print data_matrix
		nrows, ncols = data_matrix.shape
		data_sum = numpy.sum(filter_func(data_matrix), axis=0)[1:]/float(nrows)  #sum over rows of the matrix, and get rid of the first column 
		data_sum_reordered = reorder(data_sum, chain_num)
		data_all.append(data_sum_reordered)
		chain_num+=1

	return data_all

def matshow_axis():
	r=csv.reader(open('chain1_table.dat' % vars()), delimiter=' ')
	residue_map=[]
	residue_map_dict = {}
	for line in r:
		residue_map.append(line[0])
		# residue_map_dict[line[0]]=int(line[1])

	#print dictionary in sorted order
	residue_map_in_order = []
	for key in sorted(residue_map, key=lambda k: int(re.search('[0-9]+', k).group(0))):
		# print key, residue_map[key]
		residue_map_in_order.append(key[0:3])
	return residue_map_in_order


def plot_matshow(filename, data, data_axis=None):
	""" plots a matrix using the matshow function with a specific formating """

	pylab.rcParams['xtick.labelsize']='8'
	pylab.rcParams['ytick.labelsize']='8'
	pylab.rcParams['legend.fontsize']='8'	
	pylab.axis('tight')
	
	# pylab.matshow(numpy.transpose(data), cmap=pylab.cm.spectral)
	pylab.matshow(data, cmap=pylab.cm.spectral)
	pylab.clim([0, data.max()])
	pylab.colorbar(orientation="horizontal", shrink=1)
	#pylab.yticks(x, tuple(data_axis))
	
	residue_map_in_order = matshow_axis()
	print "the plot axis", residue_map_in_order
	
	x = numpy.arange(0,len(data[0]))
	pylab.xticks(x, tuple(residue_map_in_order))
	y = numpy.arange(0, len(data[:,0]))
	pylab.yticks(y, tuple(["chain 1", "chain 2", "chain 3", "chain 4", "chain 5"]))
	pylab.savefig(filename + '.pdf' % vars())	
	numpy.savetxt(filename + '.dat', numpy.transpose(data), fmt='%f')

def gen_figure(name, data):
	all_data_as_matrix = numpy.array(data)
	scyllo = numpy.average(all_data_as_matrix[0:10], axis=0)
	chiro = numpy.average(all_data_as_matrix[11:], axis=0)
	#nonpolar_axis = plot_matshow_axis()
	plot_matshow(name+'_sc', scyllo)
	plot_matshow(name+'_ch', chiro)

def main():
	"""docstring for main"""
	pass


if __name__ == '__main__':
	# system = sys.argv[1]
	h5list = glob.glob("*.h5")
	
	tables_list_nonpolar = ['/nonpolar/inositol_residue_np_ch1','/nonpolar/inositol_residue_np_ch2',
	'/nonpolar/inositol_residue_np_ch3','/nonpolar/inositol_residue_np_ch4',
	'/nonpolar/inositol_residue_np_ch5']
	
	tables_list_hbonds = ['/hbond/inositol_hbond_chain0', '/hbond/inositol_hbond_chain1', '/hbond/inositol_hbond_chain2', '/hbond/inositol_hbond_chain3', '/hbond/inositol_hbond_chain4' ]
	
	
	# print "list of h5 files to read:", h5list
	# 	print "reading tables:", tables_list
	# 	
	all_systems_nonpolar_contact = []
	all_systems_polar_contact = []
	for h5file in h5list:
		print h5file
		f = tables.openFile(h5file, mode='a')
		
		data_nonpolar = build_plot_data(f, tables_list_nonpolar)
		all_systems_nonpolar_contact.append(data_nonpolar)
		data_polar = build_plot_data(f, tables_list_hbonds, filter_func=filter_hbond)
		all_systems_polar_contact.append(data_polar)
	
	gen_figure('nonpolar', all_systems_nonpolar_contact)
	gen_figure('polar', all_systems_polar_contact)

