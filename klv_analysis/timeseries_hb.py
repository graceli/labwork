#!/usr/bin/python
import re
import glob
import os

import pylab
import numpy

import plot_and_save2hdf5 as myh5
import utils
import config

#helper functions
def create_datalist(ratio):
	"""docstring for create_datalist"""
	filenames = glob.glob("*%(ratio)s*p2p_vs_t*smoothed.txt.gz" % vars())
	assert len(filenames)>0, "expecting more than 0 files"

	datalist=[]
	labellist=[]
	if len(filenames):
		for f in filenames:
			fields = f.split('_')
			iso = fields[0]
			conc = fields[1]
			datalist.append(numpy.genfromtxt(f))
			if iso == "water":
				labellist.append(iso)
			else:
				labellist.append(iso + " " + config.RATIO[conc])

	return (datalist, labellist)

def process(h5file, ratio, format="p2p_vs_t.dat"):
	# given a h5file return a list of data to be plotted as line plots 
	# and a corresponding list of labels

	header = "# time average_inter std_inter average_intra std_intra"
	datalist = []
	labellist = []
	isomerlist = ["scyllo", "chiro", "water"]

	for iso in isomerlist:
		print "processing", iso
		pattern = re.compile(r"%(iso)s.*%(ratio)s.*%(format)s" % vars())
		if iso == "water":
			pattern = re.compile(r"%(iso)s.*%(format)s" % vars())
			
		data_inter = []
		data_intra = []
		for table in h5file.listNodes(where='/polar'):
			table_path = os.path.join('/polar', table.name)
			if pattern.search(table.name):
				print "processing", table.name
				data = myh5.getTableAsMatrix(h5file, table_path, dtype=numpy.int32)
				data = data.astype('float')
				print "converted to float32", data
				
				nrows, ncols = data.shape
				assert nrows > ncols
				print "Test data read in dimensions", data.shape, data.dtype
				data_inter.append(data[0:config.LASTFRAME,1])
				data_intra.append(data[0:config.LASTFRAME,2])
				
		# compute summary statistics
		print "summarizing statistics ... "
		
		average_inter, std_inter = utils.summary_statistics(utils.array_list_to_matrix(data_inter))
		average_intra, std_intra = utils.summary_statistics(utils.array_list_to_matrix(data_intra))
		time = data[0:config.LASTFRAME,0]
		
		# print "Test: dimensions of average_inter", average_inter.shape
		plotdata = utils.array_list_to_matrix([ time, average_inter, std_inter, average_intra, std_intra ])
		print "plotdata", plotdata
		print "Test: dimensions of plotdata for", iso, ratio, plotdata.shape
		plotdata_smoothed = utils.smooth(plotdata, 500, time_present=True, timestep=2)
		print plotdata_smoothed
		
		datalist.append(plotdata_smoothed)
		print "smoothed data", plotdata_smoothed, plotdata_smoothed.shape
		
		ratiolabel = config.RATIO[ratio]
		if iso == "water":
			labellist.append("%(iso)s" % vars())
		else:
			labellist.append("%(iso)s (%(ratiolabel)s)" % vars())
			
		utils.savetxt('%(iso)s_%(ratio)s_p2p_vs_t.txt' % vars(), header, plotdata, fmt='%0.2f')
		utils.savetxt('%(iso)s_%(ratio)s_p2p_vs_t_smoothed.txt' % vars(), header, plotdata_smoothed, fmt='%0.2f')

	return (datalist, labellist)
	
def plot(datalist, labellist, ratio, use_flat_files=False):
	if use_flat_files == True:
		print "using flat files ... "
		datalist,labellist = create_datalist(ratio)

	assert len(datalist) == 3, "%s data lists, expecting 3 for scyllo, chiro, water" % `len(datalist)`
	
	config.configure_plot()
	fig1 = pylab.figure(num=1)
	fig2 = pylab.figure(num=2)
	ax1 = fig1.add_subplot(111)
	ax2 = fig2.add_subplot(111)

	print "plotting data"
	for i in range(0,len(datalist)):
		parts = labellist[i].split()
		isomer = parts[0]
		print "plotting", isomer
		
		if isomer == "water":
			label = "no inositol"

		# print "Test: time column for", labellist[i], datalist[i][:,0]
		nrows, ncols = datalist[i].shape
		assert nrows > ncols, "Number of cols greater than rows!"

		time = datalist[i][:,0]/1000
		ax1.fill_between(time, 
						(datalist[i][:,1] - datalist[i][:,2])/config.NMOLECULES, 
						(datalist[i][:,1] + datalist[i][:,2])/config.NMOLECULES,
						alpha=1,
						facecolor=config.SHADED_COLOR[isomer],
						edgecolor=config.SHADED_COLOR[isomer], 
						lw=0.5)
		ax2.fill_between(time,
					(datalist[i][:,3] - datalist[i][:,4])/config.NMOLECULES,
					(datalist[i][:,3] + datalist[i][:,4])/config.NMOLECULES, 
					alpha=1,
					facecolor=config.SHADED_COLOR[isomer],
					edgecolor=config.SHADED_COLOR[isomer], 
					lw=0.5)
		ax1.plot(time, datalist[i][:,1]/config.NMOLECULES, color=config.LINE_COLOR[isomer], label=labellist[i])
		ax2.plot(time, datalist[i][:,3]/config.NMOLECULES, color=config.LINE_COLOR[isomer], label=labellist[i])

	ax2.set_xlabel('Time (ns)')
	ax1.set_ylabel('Intermolecular Hydrogen Bonds per peptide')
	ax2.set_ylabel('Intramolecular Hydrogen Bonds per peptide')
	ax1.grid(True); ax2.grid(True)
	ax1.set_xlim(0, config.RUNTIME_NS); ax2.set_xlim(0, config.RUNTIME_NS)
	ax1.set_ylim(0,2); ax2.set_ylim(0,2)

	ax1.legend(loc='lower right', ncol=1, columnspacing=0.5, borderaxespad=0.)
	ax2.legend(loc='lower right', ncol=1, columnspacing=0.5, borderaxespad=0.)

	#frame = leg.get_frame()
	#frame.set_linewidth(0.5)
	print "Saving figure..."

	#this produces a small image of this size, but everything else is way too big
	#need to adjust rcParams for font size and legend size
	# fig1.set_size_inches((3.25,3.5))
	# fig2.set_size_inches((3.25,3.5))
	print "Polar figure 1 size", fig1.get_size_inches()
	print "Polar figure 2 size", fig2.get_size_inches()

	fig1.savefig('inter_%(ratio)s.png' % vars())
	fig2.savefig('intra_%(ratio)s.png' % vars())


def run(h5file, ratio, use_flat_files=False):
	""" plot the peptide-peptide time series """
	# data_nonpolar = get_time_series(h5file, '/pp_nonpolar')
	# plot_timeseries(data_nonpolar)
	# 
	# data_polar = polar(h5file, '/polar')
	# plot_timeseries(data_polar)
	if use_flat_files == True:
		data_list=[]
		label_list=[]
		plot(data_list, label_list, ratio, use_flat_files)
	else:
		(data_list, label_list) = process(h5file, ratio)
		plot(data_list, label_list, ratio, use_flat_files)