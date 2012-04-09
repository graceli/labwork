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
	""" reorder the nonpolar data according to a textfile mapping """

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


def build_plot_data_nonpolar_contact(h5file, tables):
	# This converts the number of hydrogen bonds into whether it is bound or not
	# Read chain tables and reorder data according to the residue mapping
	data_all = []
	chain_num = 1
	for table in tables:
		datatable = h5file.getNode(table).read()
		row_size = len(datatable[0])
		data_matrix = datatable.view(dtype=numpy.float64).reshape(-1, row_size)

		filtered_matrix = filter_nonpolar(data_matrix)
		nrows, ncols = filtered_matrix.shape

		print " analyzing nonpolar data of dimension", filtered_matrix.shape, "from", table
		data_sum=[]
		for c in range(0, ncols):
			#calculate the number of nonzero elements in a column (a column corresponds to a residue)
			column = filtered_matrix[:,c]		
			elements = filtered_matrix[numpy.nonzero(column),c][0]
			
			#fraction bound 
			data_sum.append(len(elements)/float(nrows))
			print nrows
		
		data_sum_reordered = reorder(numpy.array(data_sum), chain_num)
		data_all.append(data_sum_reordered)
		chain_num += 1

	return data_all

def build_plot_data_polar_contact(h5file, tables):
	""" Builds a list of data for each given table to be plotted
		This changes depending on the analysis """

	#read chain tables and reorder data according to the residue mapping	
	data_all = []
	chain_num = 1
	for table in tables:
		print table
		datatable = h5file.getNode(table).read()
		row_size = len(datatable[0])
		data_matrix = datatable.view(dtype=numpy.float64).reshape(-1, row_size)
	
		filtered_matrix = filter_hbond(data_matrix)
		nrows, ncols = filtered_matrix.shape

		print " analyzing nonpolar data of dimension", filtered_matrix.shape, "from", table

		data_sum = []
		for c in range(0, ncols):
			#calculate the number of nonzero elements in a column (a column corresponds to a residue)
			column = filtered_matrix[:,c]
			elements = filtered_matrix[numpy.nonzero(column),c][0]
			print elements
			print len(elements)
			data_sum.append(len(elements)/float(nrows))
			print nrows
		
		data_all.append(numpy.array(data_sum))
		chain_num += 1

	return data_all


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

def matshow_axis_polar_full():
	""" hardcoded residue names for the polar data axis (all the residues in Abeta17-42) """

	residue_map = ['LEU', 'VAL', 'PHE', 'PHE', 'ALA', 'GLU', 'ASP', 'VAL', 'GLY', 'SER', 'ASN', 'LYS',
	'GLY', 'ALA', 'ILE', 'ILE', 'GLY', 'LEU', 'MET', 'VAL','GLY', 'GLY', 'VAL','VAL', 'ILE', 'ALA']
	
	return residue_map

def matshow_axis_polar():
	""" hardcoded residue names for the polar data axis (all the residues in Abeta17-42) """

	residue_map = ['L', 'V', 'F', 'F', 'A', 'E', 'D', 'V', 'G', 'S', 'N', 'K',
	'G', 'A', 'I', 'I', 'G', 'L', 'M', 'V','G', 'G', 'V','V', 'I', 'A']

	return residue_map
		
def plot_matshow(filename, data, data_axis=None, limits=[0,0.03]):
	""" plots a matrix using the matshow function with a specific formating """
	
	#how the matrix will be rendered
	pylab.matshow(data, cmap=pylab.cm.spectral)
	pylab.clim(limits)
	pylab.colorbar(orientation="horizontal", shrink=1)

	print "the plot axis labels:", data_axis	

	#how the axis will be rendered
	x = numpy.arange(0,len(data[0]))
	xlabels=tuple(data_axis)
	pylab.xticks(x, xlabels)

	y = numpy.arange(0, len(data[:,0]))
	ylabels=tuple(["chain 1", "chain 2", "chain 3", "chain 4", "chain 5"])
	pylab.yticks(y, ylabels)

	#save the figure to disk as a pdf
	pylab.savefig(filename + '.pdf' % vars())	
	numpy.savetxt(filename + '.dat', numpy.transpose(data), fmt='%f')


def gen_contact_figure(name, data, data_axis=None, limits=None):
	all_data_as_matrix = numpy.array(data)
	print "all_data_as_matrix dimensions", all_data_as_matrix.shape
	print all_data_as_matrix

	data = numpy.average(all_data_as_matrix, axis=0)
	# plot_matshow(name, data, data_axis=data_axis, limits=[data.min(), data.max()])
	plot_matshow(name, data, data_axis=data_axis, limits=limits)
	# plot_matshow(name+'_ch', chiro, data_axis=data_axis, limits=[data_min, data_max])

def show_contact(name, h5list):
	""" computes and plots the nonpolar/polar contacts as contact maps"""
		
	pylab.rcParams['xtick.labelsize']='9'
	pylab.rcParams['ytick.labelsize']='9'
	pylab.rcParams['legend.fontsize']='9'
	pylab.rcParams['figure.figsize'] = [4,3]	
	pylab.axis('tight')
	
	tables_list_nonpolar = ['/nonpolar/inositol_residue_np_ch1','/nonpolar/inositol_residue_np_ch2',
	'/nonpolar/inositol_residue_np_ch3','/nonpolar/inositol_residue_np_ch4', '/nonpolar/inositol_residue_np_ch5']
	
	tables_list_hbonds = ['/hbond/inositol_hbond_chain0', '/hbond/inositol_hbond_chain1', '/hbond/inositol_hbond_chain2', 
	'/hbond/inositol_hbond_chain3', '/hbond/inositol_hbond_chain4' ]
	
	all_systems_nonpolar_contact = []
	all_systems_polar_contact = []
	for h5file in h5list:
		data_nonpolar = build_plot_data_nonpolar_contact(h5file, tables_list_nonpolar)
		all_systems_nonpolar_contact.append(data_nonpolar)
			
		data_polar = build_plot_data_polar_contact(h5file, tables_list_hbonds)
		all_systems_polar_contact.append(data_polar)

	nonpolar_axis = matshow_axis_nonpolar()
	polar_axis = matshow_axis_polar()

	gen_contact_figure(name + '_nonpolar', all_systems_nonpolar_contact, data_axis=nonpolar_axis, limits=[0,0.03])
	gen_contact_figure(name + '_polar', all_systems_polar_contact, data_axis=polar_axis, limits=[0, 0.3])

def gen_rmsd_figure(nameslist, h5list):
	print len(h5list), "analyzed"
	pylab.rcParams['xtick.labelsize']='9'
	pylab.rcParams['ytick.labelsize']='9'
	pylab.rcParams['legend.fontsize']='9'
	pylab.rcParams['figure.figsize'] = [6,4]
	pylab.rcParams['figure.subplot.wspace'] = 0.25
	pylab.rcParams['figure.subplot.hspace'] = 0.3
	
	tables_list = ['/protein/rmsd_ca2']
	fig_rows=1
	fig_cols=1
	fig_num = 1
	lines_list = []
	for f in h5list:
		data = f.getNode('/protein/rmsd_ca2').read()
 		data_matrix = data.view(dtype=numpy.float64).reshape(-1, len(data[0]))
		pylab.subplot(fig_rows, fig_cols, fig_num)
		# pylab.axis('tight')
		x=data_matrix[:,0]/1000.0
		y=data_matrix[:,1]
		print x, len(x)
		print y, len(y)
		
		plist = pylab.plot(x,y, 'g', label=nameslist[fig_num-1])
		pylab.axis([0,200,0.2,0.8])
		# pylab.legend(loc=0)
		lines_list.append(plist)
		pylab.plot(x, [0.5]*len(x))
		fig_num+=1
	
	# pylab.legend(linees_list, ('%(name)s backbone' % vars(),), loc=0)
	# pylab.legend(loc=0)

	
def show_rmsd(h5list,name=''):
	# nameslist = ['scyllo', 'chiro']
	nameslist = [ 'h5' + '/' + 'water' + str(i) + '.h5' for i in range(1,11) ]
	gen_rmsd_figure(nameslist, h5list)
	# gen_rmsd_figure('chiro', h5list)
	# pylab.xlabel('Time (ns)')
	# pylab.ylabel('RMSD (nm)')
	pylab.savefig(name + '.png')
	pylab.clf()


def show_rmsf(h5list, name=''):
	pylab.rcParams['xtick.labelsize']='5'
	pylab.rcParams['ytick.labelsize']='6'
	pylab.rcParams['legend.fontsize']='6'
	pylab.rcParams['figure.figsize'] = [6,4]
	pylab.rcParams['figure.subplot.wspace'] = 0.2
	pylab.rcParams['figure.subplot.left'] = 0.08
	pylab.rcParams['figure.subplot.right'] = 0.99

	# nameslist=['scyllo', 'chiro']
	# nameslist = [ 'h5' + '/' + 'water' + str(i) + '.h5' for i in range(1,11) ]
	nameslist = [ 'glycerol' ]
	fig_rows=1
	fig_cols=1
	fig_num=1
	for f in h5list:
		table = f.getNode('/protein/rmsf_ca2')
		data = table.read().view(dtype=numpy.float64).reshape(-1, len(table[0]))
		values = data[:,1]
		
		values_chain_matrix = values.reshape(-1, len(values)/5)
		rows, cols = values_chain_matrix.shape
		# average_rmsf_per_residue = numpy.average(values_chain_matrix, axis=0)
		# pylab.axis('off')

		pylab.subplot(fig_rows, fig_cols, fig_num)
		for i in range(0, rows):
			# pylab.plot(values_chain_matrix[i], label=nameslist[fig_num-1]+' chain %(i)d' % vars())
			pylab.plot(values_chain_matrix[i], label=nameslist[fig_num-1]+' chain %(i)d' % vars())
		
		pylab.xticks(numpy.arange(0,26), matshow_axis_polar())
		pylab.yticks()
		pylab.grid(True)
		# pylab.xlabel('Abeta1-42 residues')
		# pylab.ylabel('RMSF (nm)')
		pylab.legend(loc=0)
	
		fig_num += 1
	
	pylab.axis([0,25,0,0.6])
	pylab.xticks(numpy.arange(0,26), matshow_axis_polar())
	pylab.savefig(name + '.png')

def main():
	"""docstring for main"""
	
	# #open all h5 files found
	# scyllo = openh5(['h5/scyllo.h5'])
	# chiro = openh5(['h5/chiro.h5'])
	 
	# show_contact('scyllo', scyllo)
	# show_contact('chiro', chiro)
 
	# scyllo.extend(chiro)
	fileslist = [ 'h5' + '/' + 'water' + str(i) + '.h5' for i in range(1,11) ]
	# water = openh5(fileslist)

	cer = openh5(['h5/cer.h5'])

	#show_rmsd()
	# print scyllo, len(scyllo)
	
	show_rmsd(cer, 'rmsd')
	show_rmsf(cer, 'rmsf')
	
	closeh5(cer)

	
	# closeh5(scyllo)
	# closeh5(chiro)

if __name__ == '__main__':
	main()
