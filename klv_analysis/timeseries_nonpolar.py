#!/usr/bin/env python
import numpy
import sys
import os
import re
import pylab
import config
import plot_and_save2hdf5 as myh5
import utils

def plot(data, ratio):
	#hack -- the matrix is in the format of time avg std (scyllo) time avg std (chiro) time avg std (water)
	assert data.shape[0] > data.shape[1], "nrows in data matrix should be greater than ncols"
	assert data.shape[1] == 9, "Expecting 6 columns, found only %s" % `data.shape[1]`
	
	# pylab.clf()
	print "nonpolar size", pylab.gcf().get_size_inches()

	#plot specific
	# fig = pylab.figure()
	fig = pylab.figure()
	ax = fig.add_subplot(111)
	
	ax.grid(True)
	ax.set_xlabel("Time (ns)")
	ax.set_ylabel("Intermolecular nonpolar contacts per peptide")
	ratio_label = config.RATIO[ratio]

	ax.fill_between(data[:,0]/1000, data[:,1] - data[:,2], data[:,1] + data[:,2],
	facecolor='chartreuse',
	alpha=1, edgecolor='chartreuse',
	lw=0.5, label='scyllo- err')
	ax.plot(data[:,0]/1000, data[:,1], label='scyllo %(ratio_label)s' % vars(), linewidth=1, color='green')
	
	ax.fill_between(data[:,3]/1000, data[:,4] - data[:,5], data[:,4] + data[:,5],
	facecolor='pink',
	alpha=1, edgecolor='pink',
	lw=0.5, label='chiro- err')
	ax.plot(data[:,3]/1000, data[:,4], label='chiro %(ratio_label)s' % vars(), linewidth=1, color='red')
	
	ax.fill_between(data[:,6]/1000, data[:,7] - data[:,8], data[:,7] + data[:,8],
	facecolor='lightblue',
	alpha=1, edgecolor='lightblue',
	lw=0.5, label='no inositol err')
	ax.plot(data[:,6]/1000, data[:,7], label='no inositol', linewidth=1,color='blue')

	ax.legend(loc='upper right')
	# pylab.draw()
	
	ax.set_xlim(0, config.RUNTIME_NS)
	ax.set_ylim(0, 30)
	# pylab.draw()

	print "Saving figure..."
	fig.savefig('%(ratio)s_pp_nonpolar.png' % vars())
	# pylab.clf()


def process(h5file, ratio):
	isomerlist = ["scyllo", "chiro", "water"]
	plot_data = []

	#read in files for each system and aggregate
	format="pp_nonpolar_vs_t.xvg"
	for iso in isomerlist:
		print "processing", iso
		pattern = re.compile(r"%(iso)s.*%(ratio)s.*%(format)s" % vars())
		if iso == "water":
			pattern = re.compile(r"%(iso)s.*%(format)s" % vars())

		datalist=[]
		for table in h5file.listNodes(where='/pp_nonpolar'):
			table_path = os.path.join('/pp_nonpolar', table.name)
			if pattern.search(table.name):			
				data = myh5.getTableAsMatrix(h5file, table_path)
				if data is not None:
					data = data.astype('float')
					datalist.append(data[0:config.LASTFRAME, 1])
				else:
					print "no data was read in"
			
		print "datalist", datalist
		data_matrix = numpy.transpose(numpy.vstack(datalist))	
		print "data_matrix", data_matrix, data_matrix.shape

		avg, std = utils.summary_statistics(data_matrix, sum_across="columns")
		avg_smoothed = utils.smooth(avg/config.NMOLECULES, 500, time_present=False, timestep=2)
		std_smoothed = utils.smooth(std/config.NMOLECULES, 500, time_present=True, timestep=2)
		plot_data.append(avg_smoothed)
		plot_data.append(std_smoothed)
	
	timeseries_matrix = numpy.hstack(plot_data)
	print "timeseries_matrix", timeseries_matrix, timeseries_matrix.shape
	print "time", timeseries_matrix[:,0]
	numpy.savetxt(ratio + "_pp_nonpolar_smoothed.txt.gz", timeseries_matrix, fmt='%0.3f')
	return timeseries_matrix

def run(h5file, ratio, use_flat_files=False):
	if use_flat_files:
		filename = ratio + "_pp_nonpolar_smoothed.txt.gz"
		data = numpy.genfromtxt(filename, comments="#")
		assert data.shape[0] > data.shape[1]
		plot(data, ratio)
	else:
		data = process(h5file, ratio)
		assert data.shape[0] > data.shape[1]
		plot(data, ratio)
