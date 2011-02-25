import tables
import numpy
import pylab
import os
# import glob
import re

import plot_and_save2hdf5 as myh5
import utils
import config

# def load(filename, tablename):
# 	# table_descr = {'time':tables.Float64Col(pos=0), 'clustsize':tables.Float64Col(pos=1)}
# 	file = tables.openFile('clustsize.h5',mode='a', filters=tables.Filters(complevel=8, complib='zlib'))
# 	
# 	if not file.__contains__('/clustsize_array'):
# 		group = file.createGroup(file.root, 'clustsize_array', filters=tables.Filters(complevel=8, complib='zlib'))
# 	else:
# 		group = file.getNode('/clustsize_array')
# 		
# 	data = numpy.genfromtxt(filename, comments="#")
# 	#Note that arrays don't support compression ! but savings come from binary storage on disk
# 	# but for my data files ... its the same size and a bit more on disk as hdf5
# 	file.createArray(group, tablename, data)
# 	file.flush()
# 	file.close()
		
# def plot_from_h5():
# 	if os.path.exists('clustsize.h5'):
# 		file = tables.openFile('clustsize.h5', mode='a')
# 		print "loading data"
# 		data_scyllo_15to4 = numpy.array(file.root.clustsize_array.scyllo_15to4)
# 		data_chiro_15to4 = numpy.array(file.root.clustsize_array.chiro_15to4)
# 		data_scyllo_45to4 = numpy.array(file.root.clustsize_array.scyllo_45to4)
# 		data_chiro_45to4 = numpy.array(file.root.clustsize_array.chiro_45to4)
# 
# 		# data_chiro_15to4 = file.root.clustsize.chiro_15to4.read()
# 		# data_scyllo_45to4 = file.root.clustsize.scyllo_45to4.read()
# 		# data_chiro_45to4 = file.root.clustsize.chiro_45to4.read()
# 		# pylab.plot(data_scyllo_15to4, data_chiro_15to4, data_scyllo_45to4, data_chiro_45to4)
# 		print "plotting"
# 		pylab.rcParams['font.size']=9
# 		pylab.rcParams['font.family']='Arial'
# 		fig = pylab.figure(num=None, figsize=(8, 6))
# 		ax = fig.add_subplot(111)
# 		ax.plot(data_scyllo_15to4[:,0]/1000.0, data_scyllo_15to4[:,1], label='scyllo (15 to 4)')
# 		ax.plot(data_chiro_15to4[:,0]/1000.0, data_chiro_15to4[:,1], label='chiro (15 to 4)')
# 		ax.plot(data_scyllo_45to4[:,0]/1000.0, data_scyllo_45to4[:,1], label='scyllo (45 to 4)')
# 		ax.plot(data_chiro_45to4[:,0]/1000.0, data_chiro_45to4[:,1], label='chiro(45 to 4)')
# 
# 		ax.set_xlabel('Time (ps)')
# 		ax.set_ylabel('Number of peptide clusters')
# 		ax.grid(True)
# 		#ax.set_title('Number of clusters vs time')
# 		ax.legend()
# 		print "saving fig..."
# 
# 		# Figure.set_figsize_inches( (w,h) )
# 		fig.savefig('cluster_large.png')
# 		#this produces a small image of this size, but everything else is way too big
# 		#need to adjust rcParams for font size and legend size
# 
# 		fig.set_size_inches((3.25,3.25))
# 		fig.savefig('cluster_small.png')

def plot(name, datalist, labellist):
		print "plotting from array"
		print len(datalist)
		print labellist
		
		pylab.rcParams['font.size']=8
		pylab.rcParams['font.family']='Arial'
		pylab.rcParams['legend.fontsize']='small'
		
		fig = pylab.figure(num=None)
		fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95, wspace=0, hspace=0)
		
		ax = fig.add_subplot(111)
		
		for i in range(0,len(datalist)):
			ax.plot(datalist[i][:,0]/1000.0, datalist[i][:,1], label=labellist[i], color=config.LINE_COLOR[config.ISOMER_LIST[i]])

		ax.set_xlabel('Time (ns)')
		ax.set_ylabel('Number of peptide clusters')
		ax.grid(True)
		ax.set_xlim(0, config.RUNTIME_NS)
		ax.legend(loc='upper right')
		print "saving fig..."
	
		#this produces a small image of this size, but everything else is way too big
		#need to adjust rcParams for font size and legend size
		fig.set_size_inches((3.25,3.25))
		pylab.draw()
		fig.savefig(name)

def process(h5file, ratio):
	labellist = []
	plot_list = []
	isomerlist = [ "scyllo", "chiro", "water" ]
	root='/cluster'
	pattern=""	
	for iso in isomerlist:
		datalist=[]
		for t in h5file.listNodes(where=root):
			filename = t.name
			# fileslist = glob.glob("*%(iso)s*%(conc)s*.xvg" % vars())
			if iso != "water":
				pattern = re.compile(r'%(iso)s.*%(ratio)s.*nclust.xvg' % vars())
			else:
				pattern = re.compile(r'%(iso)s.*nclust.xvg' % vars())
		
			if pattern.search(filename):
				print filename
				data = myh5.getTableAsMatrix(h5file, 
									 os.path.join(root, filename))[0:config.LASTFRAME]
				print data.shape
				datalist.append(data)
		
		print len(datalist)
		all_data = numpy.hstack(datalist)
		print all_data.shape
		print all_data
			
		sdata = all_data[:,1::2]			
		average = numpy.average(sdata, axis=1)
		plot_data = numpy.transpose([all_data[:,0],average])
		
		numpy.savetxt('%(iso)s_%(ratio)s_nclust.txt.gz' % vars(), plot_data, fmt='%0.2f')
		
		plot_list.append(plot_data)
		plot_label = config.LABEL[iso] + " (" + config.RATIO[ratio] + ")"
		
		labellist.append("%(plot_label)s" % vars())
		 
	return (plot_list, labellist)
				
#load('chiro_130ns_15to4_nclust.txt', 'chiro_15to4')
#load('scyllo_130ns_15to4_nclust.txt', 'scyllo_15to4')
#load('water_130ns_nclust.txt', 'water')

# (datalist, labellist) = process()
# plot(datalist, labellist)

def run(h5file, ratio, use_flat_files=False):
	(plot_list, labellist) = process(h5file, ratio)
	plot(ratio + '_cluster' + '.pdf', plot_list, labellist)
