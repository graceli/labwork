#!/usr/bin/env python
import pylab
import numpy
import matplotlib
import sys

def plot_settings(setting="default"):
	if setting == "pub":
		pylab.rcParams['figure.figsize'] = [6,3]
		pylab.rcParams['font.size'] = 9
		pylab.rcParams['legend.fontsize'] = 'small'
	else:
		pylab.rcParams['figure.figsize'] = [8,6]
		pylab.rcParams['font.size'] = 12


def plot_with_error(x, y, yerror, label, axes=None):
	""" plots a figure onto a figure object and 
	returns the modified figure"""
	
	axes.fill_between(x, y-yerror, y+yerror, alpha=0.1)
	axes.plot(x,y,label=label)
	# fig.savefig(name)
	
	return axes

def plot_eed(file_name="test.png"):
	# I don't have Latex installed
	# pylab.rcParams['text.usetex']=True
	# pylab.rcParams['text.latex.unicode']=True
	
	scyllo = numpy.genfromtxt('scyllo_eed_data2.dat', dtype=float)	
	chiro = numpy.genfromtxt('chiro_eed_data2.dat', dtype=float)	
	plot_settings(setting="pub")
	
	fig = pylab.figure()
	# fig.subplots_adjust(top=0.9, bottom=0.19, left=0.18, right=0.90, wspace=0.001)
	pylab.axes([0.05, 0.05, 0.92, 0.92])

	ax1 = fig.add_subplot(1,2,1)
	# figure sharing doesn't work -- don't know why ...
	ax2 = fig.add_subplot(1,2,2)	
	ax1.set_xlabel(r'End to end distance ($\AA$)')
	ax1.set_ylabel("Probability")
	yticklabels = ax2.get_yticklabels()+ax2.get_xticklabels()
	pylab.setp(yticklabels, visible=False)

	nax = 0
	for eed_data in [scyllo, chiro]:
		nrows, ncols = eed_data.shape
		ax = fig.get_axes()[nax]
		ax.grid('on')
		ax.set_xlim(0, 25)
		ax.set_ylim(0, 0.1)
		for i in range(1, ncols, 2):
			print 0, i, i+1
			x = eed_data[:,0]
			y = eed_data[:,i]
			yerror = eed_data[:,i+1]
			state = int(i/2)
			ax = plot_with_error(x, y, yerror, '%(state)d hbonds' % vars(), axes=ax)
		nax += 1
	ax2.legend()
	fig.savefig("all.png")

def plot_timeseries(xlabel='', ylabel=''):
	""" plot the time series for peptide-peptide nonpolar and polar contacts """
	
	pylab.rcParams['figure.figsize'] = [2.9, 2.75]
	pylab.rcParams['font.size'] = 9
	pylab.rcParams['legend.fontsize'] = 'small'
	pylab.rcParams['legend.loc'] = 'best'
	# pylab.subplots_adjust(top=0.9, bottom=0.19)
	# pylab.axes([0.05, 0.05, 0.92, 0.92])
	
	for system in ['systematic_4oct', 'systematic_ap1f', 'systematic_ap2f']:
		for iso in ['scyllo', 'chiro', 'control']:
			filename = '%(system)s_%(iso)s_average.txt' % vars()
			try:
				average = numpy.genfromtxt(filename)
			except IOError:
				print filename,"is empty"
			else:
				pylab.axes([0.15, 0.12, 0.81, 0.85])
				pylab.xlim(0, 100)
				pylab.ylim(0, 13)
				
				x = average[:,0]/1000
				y1 = average[:,1]
				y2 = average[:,2]
			
				pylab.plot(x, y1, label='%(iso)s' % vars(), linewidth=1)
				pylab.plot(x, y2, label='%(iso)s' % vars(), linewidth=1)
			
				pylab.grid(True)
				# pylab.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper right',
				       # ncol=3, borderaxespad=0.)
				pylab.legend(loc='upper center', ncol=3, borderaxespad=0)


		pylab.ylabel(ylabel)
		pylab.xlabel(xlabel)
		pylab.savefig('%(system)s.png' %vars())
		pylab.clf()

def plot_timeseries_single(xlabel='', ylabel=''):
	""" plot the time series for peptide-peptide nonpolar and polar contacts """

	pylab.rcParams['figure.figsize'] = [2.75,2.75]
	pylab.rcParams['font.size'] = 9
	pylab.rcParams['legend.fontsize'] = 'small'
	pylab.rcParams['legend.loc'] = 'best'

	for system in ['systematic_4oct', 'systematic_ap1f', 'systematic_ap2f']:
		for iso in ['scyllo', 'chiro', 'control']:
			filename = '%(system)s_%(iso)s_average.txt'
			try:
				average = numpy.genfromtxt('%(system)s_%(iso)s_average.txt' % vars())
			except IOError:
				print filename,"is empty"
			else:
				pylab.axes([0.12, 0.12, 0.85, 0.86])
				x = average[:,0]/1000
				y1 = average[:,1]
				# y2 = average[:,2]
				pylab.xlim(0, 100)
				pylab.ylim(0, 7)
				
				pylab.plot(x, y1, label='%(iso)s' % vars(), linewidth=1)
				# pylab.plot(x, y2, label='%(iso)s inter-' % vars(), linewidth=1)

				pylab.grid(True)
				# pylab.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
				       # ncol=3, mode="expand", borderaxespad=0.)
				pylab.legend(loc='lower right')
		
		pylab.ylabel(ylabel)
		pylab.xlabel(xlabel)
		pylab.savefig('%(system)s.png' %vars())
		pylab.clf()

# Decided not to use this to plot stuff
# def plot_vertical_bars(name, data_columns, data_error, ylabel='y-axis', xlabel='x-axis', xticks=None, legend=None):
# 	"""
# 		data_columns -- data to plot as bars organized in columns, each column is a set of bars, with the height of the bar
# 		in each row of the column.  The number of columns is the number of different category of bars.
# 	"""
# 	# TODO: assert that the data_column is a numpy array
# 	barwidth = 0.35
# 	nrows, ncols = data_columns.shape
# 	# indices used for the ticks of each bar plot
# 	index = numpy.arange(nrows)
# 	for i in data_columns:
# 		rects = pylab.bar(index, data_columns[i], barwidth, yerr=data_error[i])
# 		rectList.append(rects)
# 	
# 	#set plot legend and labels, if any
# 	pylab.ylabel(ylabel)
# 	pylab.xlabel(xlabel)
# 	pylab.xticks(index+barwidth, xticks)
# 	pylab.legend(tuple(rectList), legend)
# 	pylab.savefig(name + '.png')
# 	pylab.clf()
	
def main():
	# plot the eed data for monomers
	# plot_eed()
	
	# plot polar contact
	if sys.argv[1] == "polar":
		print sys.argv[0]
		plot_timeseries(xlabel='Time (ns)', ylabel='Total polar contact')
	
	# plot nonpolar contact plots
	if sys.argv[1] == "nonpolar":
		plot_timeseries_single(ylabel='Total nonpolar contact', xlabel='Time (ns)')

if __name__ == '__main__':
	main()
