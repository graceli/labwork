#!/usr/bin/env python
import pylab
import numpy
import csv
import re
import glob
import tables
import sys

# Not a well designed plot script

def openh5(h5list):
	openFiles = []
	for h5file in h5list:
		openFiles.append(tables.openFile(h5file, mode='a'))
	return openFiles

def closeh5(h5_reflist):
	for h5file in h5_reflist:
		h5file.close()

def filter_nonpolar(data_matrix):
	return data_matrix[:,1:]

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
		nonpolar_dict[residue_map[i]] = data[:,i]

	data_in_order = []
	for key in sorted(nonpolar_dict.keys(),key=lambda k: int(re.search('[0-9]+', k).group(0))):
		#return data in the order of residue number and normalize by the number of atoms 
		data_in_order.append(nonpolar_dict[key]/float(residue_map_dict[key]))

	return data_in_order

def matshow_axis_nonpolar():
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

def matshow_axis_polar():
	""" hardcoded residue names for the polar data axis (all the residues in Abeta17-42) """

	residue_map = ['LEU', 'VAL', 'PHE', 'PHE', 'ALA', 'GLU', 'ASP', 'VAL', 'GLY', 'SER', 'ASN', 'LYS',
	'GLY', 'ALA', 'ILE', 'ILE', 'GLY', 'LEU', 'MET', 'VAL','GLY', 'GLY', 'VAL','VAL', 'ILE', 'ALA']

	return residue_map

def calc_bound(data_all, region, chains):
	# data_all is a matrix of dimension 5 by 21 by number of snapshots
	# b1 
	beta1_bound = numpy.sum(data_all[:,region,:], axis=1)
	print "beta1_bound", beta1_bound.shape
	
	if len(chains) == 1:
		return beta1_bound[chains[0],:]
	
	beta1_bound_allch = numpy.sum(beta1_bound[chains,:], axis=0)
	print beta1_bound_allch.shape	
	
	#beta1_bound_allch should be a 2D vector with rows 1s and 0s and columns as time
	return beta1_bound_allch
	
def calc_fraction_unbound(data):
	"""data is a vector of bound and unbound states v t"""
	nonzero_indices = numpy.nonzero(data)
	num_bound = len(data[nonzero_indices])
	total = len(data)
	fraction  =  float(total-num_bound)/ num_bound
	return [num_bound, total-num_bound, fraction, fraction*49]

def build_plot_data_contact(h5file, tables_nonpolar, tables_polar):
	""" Builds a list of data for each given table to be plotted
		This changes depending on the analysis """

	#read chain tables and reorder data according to the residue mapping	
	data_all = []
	chain_num=1
	for table in tables_nonpolar:
		print table
		datatable = h5file.getNode(table).read()
		# print datatable.dtype
		row_size = len(datatable[0])
		# print "build_plot_data:", row_size
		data_matrix = datatable.view(dtype=numpy.float64).reshape(-1, row_size)
		# print data_matrix

		filtered_matrix = filter_nonpolar(data_matrix)
		nrows, ncols = filtered_matrix.shape
		data_sum=[]
		
		
		reordered_matrix = reorder(filtered_matrix, chain_num)
		data_all.append(reordered_matrix)
		chain_num+=1
	
	data_all = numpy.array(data_all)
	beta1 = calc_bound(data_all, range(1,8), range(0,5))
	turn = calc_bound(data_all, range(9,11), range(0,5))
	beta2 = calc_bound(data_all, range(12,19), range(0,5))
	NC = calc_bound(data_all, tuple((0,20)), range(0,5)) 
	edge1 = calc_bound(data_all, range(0,20), [0]) 
	edge2 = calc_bound(data_all, range(0,20), [4])
	
	
	
	data_all_polar = []
	for table in tables_polar:
		print table
		datatable = h5file.getNode(table).read()
		# print datatable.dtype
		row_size = len(datatable[0])
		# print "build_plot_data:", row_size
		data_matrix = datatable.view(dtype=numpy.float64).reshape(-1, row_size)
		# print data_matrix

		filtered_matrix = filter_hbond(data_matrix)
		nrows, ncols = filtered_matrix.shape
		# reordered_matrix = reorder(filtered_matrix, chain_num)
		data_all_polar.append(numpy.transpose(filtered_matrix))
	
	data_all_polar = numpy.array(data_all_polar)
	
	beta1_polar = calc_bound(data_all_polar, range(1,9), range(0,5))
	turn_polar = calc_bound(data_all_polar, range(10,13), range(0,5))
	beta2_polar = calc_bound(data_all_polar, range(14,24), range(0,5))
	NC_polar = calc_bound(data_all_polar, tuple((0,25)), range(0,5))
	edge1_polar = calc_bound(data_all_polar, range(0,25), [0]) 
	edge2_polar = calc_bound(data_all_polar, range(0,25), [4])
	
	print beta1.shape, beta1_polar.shape
	beta1_sum = beta1 + beta1_polar
	turn_sum = turn + turn_polar 
	beta2_sum = beta2 + beta2_polar
	NC_sum = NC + NC_polar
	edge1_sum = edge1 + edge1_polar
	edge2_sum = edge2 + edge2_polar

	
	print calc_fraction_unbound(beta1_sum)
	print calc_fraction_unbound(turn_sum)
	print calc_fraction_unbound(beta2_sum)
	print calc_fraction_unbound(NC_sum)
	print calc_fraction_unbound(edge1_sum)
	print calc_fraction_unbound(edge2_sum)
	
	all_nonpolar = numpy.sum(numpy.sum(data_all,axis=0),axis=0)
	all_polar = numpy.sum(numpy.sum(data_all_polar,axis=0),axis=0)
	
	print calc_fraction_unbound(all_nonpolar+all_polar)
	
def main():
	"""docstring for main"""

	scyllo = openh5(['scyllo.h5'])
	chiro = openh5(['chiro.h5'])

	tables_list_nonpolar = ['/nonpolar/inositol_residue_np_ch1','/nonpolar/inositol_residue_np_ch2',
	'/nonpolar/inositol_residue_np_ch3','/nonpolar/inositol_residue_np_ch4', '/nonpolar/inositol_residue_np_ch5']
	
	tables_list_hbonds = ['/hbond/inositol_hbond_chain0', '/hbond/inositol_hbond_chain1', '/hbond/inositol_hbond_chain2', 
	'/hbond/inositol_hbond_chain3', '/hbond/inositol_hbond_chain4' ]
	
	build_plot_data_contact(scyllo[0], tables_list_nonpolar, tables_list_hbonds)
	build_plot_data_contact(chiro[0], tables_list_nonpolar, tables_list_hbonds)

	closeh5(scyllo)
	closeh5(chiro)

if __name__ == '__main__':
	main()
